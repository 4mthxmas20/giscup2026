// Headless harness that runs the OFFICIAL GIS Cup 2026 evaluator against our
// submission, so we can confirm every claimed building survives their geometry
// engine (ArcGIS relate, 1 mm spatial tolerance) rather than only our own.
//
// Copied into the evaluator checkout's benchmarks/ directory by
// scripts/run_official_eval.sh, which supplies DATASET and SUBMISSION.
import { readFileSync, writeFileSync } from "node:fs";
import { resolve } from "node:path";
import { performance } from "node:perf_hooks";

import { expect, test } from "vitest";

import { parseBuildingDatasetText } from "../src/core/dataset-loader.js";
import { evaluateValidatedSubproblemAsync } from "../src/core/evaluation-engine.js";
import { parseSolutionText } from "../src/core/solution-parser.js";
import { validateSubproblemInput } from "../src/core/submission-validator.js";

test("official evaluator agrees with our claimed service scores", async () => {
  const datasetPath = resolve(process.env.DATASET ?? "datasets/GIS-cup-sample-dataset.geojson");
  const submissionPath = resolve(process.env.SUBMISSION ?? "solutions.txt");
  const reportPath = resolve(process.env.REPORT ?? "official-eval-report.json");

  const dataset = parseBuildingDatasetText(readFileSync(datasetPath, "utf8"));
  const parsed = parseSolutionText(readFileSync(submissionPath, "utf8"));
  // Optional 1-based filter so CI can fan the nine sub-problems out in parallel.
  const only = process.env.SUBPROBLEM === undefined
    ? undefined
    : Number(process.env.SUBPROBLEM);
  const selected = only === undefined
    ? parsed.subproblems
    : parsed.subproblems.filter((subproblem) => subproblem.index === only);
  if (selected.length === 0) throw new Error(`No sub-problem matched SUBPROBLEM=${only}`);
  console.log(
    `\nDATASET ${datasetPath}\n  buildings=${dataset.buildings.length} edges=${dataset.edgeCount}`
    + `\nSUBMISSION ${submissionPath}\n  subproblems=${parsed.subproblems.length}`
    + ` evaluating=${selected.length}`,
  );

  const rows: unknown[] = [];
  let totalClaimed = 0;
  let totalVerified = 0;

  for (const source of selected) {
    const validated = validateSubproblemInput(source, dataset);
    const started = performance.now();
    const result = await evaluateValidatedSubproblemAsync(dataset, validated, {});
    const seconds = (performance.now() - started) / 1000;

    const claimed = validated.claims.reportedIds.length;
    const verified = result.verifiedServiceScore;
    totalClaimed += claimed;
    totalVerified += verified;

    const warningCounts: Record<string, number> = {};
    for (const warning of result.warnings) {
      warningCounts[warning.code] = (warningCounts[warning.code] ?? 0) + 1;
    }
    const rejected = result.warnings
      .filter((warning) => warning.code === "CLAIM_BELOW_THRESHOLD")
      .slice(0, 5)
      .map((warning) => warning.message);

    const row = {
      tau: source.tau,
      k: source.k,
      antennasParsed: source.antennas.length,
      antennasRetained: source.retainedAntennas.length,
      antennasValid: validated.validAntennas.length,
      antennasUnique: validated.uniqueAntennas.length,
      antennasInvalid: validated.invalidRetainedAntennaCount,
      claimed,
      verified,
      lost: claimed - verified,
      seconds: Math.round(seconds * 10) / 10,
      warningCounts,
      sampleRejections: rejected,
    };
    rows.push(row);
    console.log(`RESULT ${JSON.stringify(row)}`);
  }

  writeFileSync(reportPath, JSON.stringify({ datasetPath, submissionPath, rows }, null, 2));
  console.log(
    `\nTOTAL claimed=${totalClaimed} verified=${totalVerified} lost=${totalClaimed - totalVerified}`
    + `\nreport written to ${reportPath}`,
  );

  // Any lost claim means we are over-reporting and would be penalized.
  expect(totalVerified).toBe(totalClaimed);
});
