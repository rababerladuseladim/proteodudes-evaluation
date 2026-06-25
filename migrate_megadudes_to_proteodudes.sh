#!/usr/bin/env bash
# Migrate on-disk workflow outputs from megadudes → proteodudes naming.
# Run from the repo root. Safe to re-run: skips moves where the source does not exist.
set -euo pipefail

move_if_exists() {
    local src="$1" dst="$2"
    [ -e "$src" ] || return 0
    mkdir -p "$(dirname "$dst")"
    mv "$src" "$dst"
    echo "Moved: $src → $dst"
}

rename_glob_in_dir() {
    # rename_glob_in_dir <dir> <old_prefix> <new_prefix>
    local dir="$1" old="$2" new="$3"
    [ -d "$dir" ] || return 0
    find "$dir" -maxdepth 1 -name "${old}*" | while read -r f; do
        base=$(basename "$f")
        new_base="${new}${base#"$old"}"
        [ "$base" != "$new_base" ] || continue
        mv "$f" "$dir/$new_base"
        echo "Renamed: $f → $dir/$new_base"
    done
}

# 1. QC plot log: logs/megadudes/plot_megadudes_qc-*.txt
#    → logs/qc/taxonomic_profiling/plot_qc_taxonomic_profiling-*.txt
for f in logs/megadudes/plot_megadudes_qc-*.txt; do
    [ -e "$f" ] || break
    base=$(basename "$f")
    move_if_exists "$f" "logs/qc/taxonomic_profiling/${base/plot_megadudes_qc-/plot_qc_taxonomic_profiling-}"
done

# 2. QC SVG: plots/megadudes/qc-*.svg → plots/taxonomic_profiling/qc-*.svg
for f in plots/megadudes/qc-*.svg; do
    [ -e "$f" ] || break
    move_if_exists "$f" "plots/taxonomic_profiling/$(basename "$f")"
done

# 3. Remaining tool-score plots: plots/megadudes/ → plots/proteodudes/
if [ -d "plots/megadudes" ]; then
    if [ -d "plots/proteodudes" ]; then
        echo "Warning: plots/proteodudes already exists; moving contents"
        find plots/megadudes -mindepth 1 -maxdepth 1 | while read -r item; do
            mv "$item" "plots/proteodudes/"
            echo "Moved: $item → plots/proteodudes/"
        done
        rmdir plots/megadudes 2>/dev/null && echo "Removed empty dir: plots/megadudes"
    else
        mv plots/megadudes plots/proteodudes
        echo "Renamed: plots/megadudes → plots/proteodudes"
    fi
fi

# 4. Results: results/megadudes/ → results/proteodudes/
if [ -d "results/megadudes" ]; then
    if [ -d "results/proteodudes" ]; then
        echo "Warning: results/proteodudes already exists; moving contents"
        find results/megadudes -mindepth 1 -maxdepth 1 | while read -r item; do
            mv "$item" "results/proteodudes/"
            echo "Moved: $item → results/proteodudes/"
        done
        rmdir results/megadudes 2>/dev/null && echo "Removed empty dir: results/megadudes"
    else
        mv results/megadudes results/proteodudes
        echo "Renamed: results/megadudes → results/proteodudes"
    fi
fi

# 5. Logs dir: logs/megadudes/ → logs/proteodudes/ (QC log already moved in step 1)
if [ -d "logs/megadudes" ]; then
    if [ -d "logs/proteodudes" ]; then
        echo "Warning: logs/proteodudes already exists; moving contents"
        find logs/megadudes -mindepth 1 -maxdepth 1 | while read -r item; do
            mv "$item" "logs/proteodudes/"
            echo "Moved: $item → logs/proteodudes/"
        done
        rmdir logs/megadudes 2>/dev/null && echo "Removed empty dir: logs/megadudes"
    else
        mv logs/megadudes logs/proteodudes
        echo "Renamed: logs/megadudes → logs/proteodudes"
    fi
fi

# 6. Rename log files within logs/proteodudes/
rename_glob_in_dir "logs/proteodudes" "build_megadudes_db-" "build_proteodudes_db-"
rename_glob_in_dir "logs/proteodudes" "run_megadudes-" "run_proteodudes-"

# 7. Benchmarks (may not exist yet)
rename_glob_in_dir "benchmarks" "build_megadudes_db-" "build_proteodudes_db-"
rename_glob_in_dir "benchmarks" "run_megadudes-" "run_proteodudes-"

# 8. QC calculation logs: logs/qc/calculate_megadudes_qc_metrics_* (may not exist yet)
for f in logs/qc/calculate_megadudes_qc_metrics_*; do
    [ -e "$f" ] || break
    base=$(basename "$f")
    move_if_exists "$f" "logs/qc/${base/calculate_megadudes_qc_metrics_/calculate_qc_metrics_for_taxonomic_profiling_}"
done

echo ""
echo "Migration complete."
echo "Next: run 'snakemake --touch' to update Snakemake's modification-time tracking,"
echo "then verify with 'snakemake -n' that no rules are unexpectedly scheduled."
