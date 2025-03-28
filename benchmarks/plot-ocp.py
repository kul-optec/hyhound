import contextlib
import json
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

filename = "ocp.json"
with contextlib.suppress(IndexError):
    filename = sys.argv[1]

# Load the JSON data
with open(filename) as f:
    data = json.load(f)

stat = "median"
metric = "real_time"  # or "cpu_time"

# Create a DataFrame from the benchmarks
df = pd.DataFrame(data["benchmarks"])
df = df[df["aggregate_name"].isna()]  # Drop aggregates generated by gbench
df_runs = df[["run_name", "real_time", "cpu_time"]].groupby("run_name")
df = df_runs.aggregate([stat, "min", "max"])

# Extract the ny parameter for the x-axis
df["run_name"] = df.index
df["ny"] = df["run_name"].str.extract(r"/ny:(\d+)(?:/|$)").astype(int)
# Extract the name
df["func_name"] = df["run_name"].apply(lambda x: x.split("/", 1)[0])
df = df.sort_values(by=["func_name", "ny"])
df["cpu_usage"] = df["cpu_time"][stat] / df["real_time"][stat]
del df["run_name"]
print(df)
functions = df["func_name"].unique()
# Prepare data for plotting
nys = df["ny"].drop_duplicates()
ny = np.unique(nys)
unit = data["benchmarks"][0]["time_unit"]

# Plotting options
plot_opts = {
    "bm_factor_schur": dict(
        label="Factor (Schur)",
        linestyle="--",
        color="tab:blue",
        marker="x",
        markersize=4,
        alpha=0.8,
        linewidth=0.8,
    ),
    "bm_solve_schur": dict(
        label="Solve (Schur)",
        linestyle="--",
        color="tab:orange",
        marker="none",
        alpha=0.8,
        linewidth=0.8,
    ),
    "bm_update_schur": dict(
        label="Update (Schur)",
        linestyle="--",
        color="tab:green",
        marker="x",
        markersize=4,
        alpha=0.8,
        linewidth=0.8,
    ),
    "bm_factor_riccati": dict(
        label="Factor (Riccati)",
        linestyle="-",
        color="tab:blue",
        marker=".",
        markersize=7,
        mfc="white",
    ),
    "bm_solve_riccati": dict(
        label="Solve (Riccati)",
        linestyle="-",
        color="tab:orange",
        marker="none",
    ),
    "bm_update_riccati": dict(
        label="Update (Riccati)",
        linestyle="-",
        color="tab:green",
        marker=".",
        markersize=7,
        mfc="white",
    ),
}

fig, ax = plt.subplots(1, 1)
for function, opts in plot_opts.items():
    function_df = df[df["func_name"] == function]
    if function_df.empty:
        continue
    print(function)
    ax.plot(function_df["ny"], 1e-3 * function_df[metric][stat].array, **opts)
    ax.fill_between(
        function_df["ny"],
        1e-3 * function_df[metric]["min"].array,
        1e-3 * function_df[metric]["max"].array,
        **{k: v for k, v in opts.items() if k in {"color"}},
        alpha=0.25
    )
ax.legend(loc="upper left")
ax.set_title("OCP factorization, solution and update run times")
ax.set_xlabel("Number of inequality constraints $n_c$")
ax.set_ylabel(r"Run time [$\mu\mathrm{s}$]")
assert unit == "ns"
ax.set_ylim(0, None)
plt.savefig(filename + ".timings.pdf")

plt.show()
