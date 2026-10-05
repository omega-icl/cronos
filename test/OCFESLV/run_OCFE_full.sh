#!/usr/bin/env bash

# Build and run every OCFE_* target defined in the Makefile.
#
# Detailed compilation and execution output is written only to a unique log.
# The terminal displays one concise status line per test case.

set -u

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$script_dir" || exit 1

sweep_log=$(mktemp \
    "OCFE_full_sweep_$(date '+%Y%m%d_%H%M%S')_XXXXXX.log")

# Determine the Makefile name.
if [[ -f Makefile ]]; then
    makefile=Makefile
elif [[ -f makefile ]]; then
    makefile=makefile
else
    echo "Error: no Makefile or makefile found in $script_dir" >&2
    exit 1
fi

# Extract explicitly defined OCFE_* targets.
mapfile -t executables < <(
    awk -F: '
        /^OCFE_.*[[:space:]]*:/ {
            target = $1
            sub(/[[:space:]]+$/, "", target)

            if (!seen[target]++)
                print target
        }
    ' "$makefile"
)

if ((${#executables[@]} == 0)); then
    echo "Error: no OCFE_* targets found in $makefile" >&2
    exit 1
fi

# Write the sweep header to the log only.
{
    echo "============================================================"
    echo "OCFE PDE sweep"
    echo "Started:   $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "Directory: $script_dir"
    echo "Makefile:  $makefile"
    echo "Targets:   ${#executables[@]}"
    echo "============================================================"
    printf '  %s\n' "${executables[@]}"
    echo
} >> "$sweep_log"

#rm *.o
echo "Log file: $sweep_log"
printf "Compiling %d targets ... " "${#executables[@]}"

# Send all make output to the sweep log only.
{
    echo "============================================================"
    echo "Compilation"
    echo "Started: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "============================================================"
} >> "$sweep_log"

make -f "$makefile" -j"$(nproc)" "${executables[@]}" \
    >> "$sweep_log" 2>&1

make_status=$?

{
    echo
    echo "Compilation finished: $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "Compilation exit status: $make_status"
    echo
} >> "$sweep_log"

if ((make_status != 0)); then
    echo "FAILED"
    echo "Compilation failed with exit status $make_status." >&2
    echo "See $sweep_log for details." >&2
    exit "$make_status"
fi

echo "OK"

failed_executables=()

for executable in "${executables[@]}"; do
    start_seconds=$SECONDS

    {
        echo "============================================================"
        echo "Running: $executable"
        echo "Started: $(date '+%Y-%m-%d %H:%M:%S %Z')"
        echo "============================================================"
    } >> "$sweep_log"

    if [[ ! -x "$executable" ]]; then
        status=126

        echo "Error: ./$executable does not exist or is not executable." \
            >> "$sweep_log"
    else
        "./$executable" >> "$sweep_log" 2>&1
        status=$?
    fi

    elapsed=$((SECONDS - start_seconds))

    {
        echo
        echo "Finished: $(date '+%Y-%m-%d %H:%M:%S %Z')"
        echo "Elapsed: ${elapsed} s"
        echo "Exit status: $status"
        echo
    } >> "$sweep_log"

    if ((status == 0)); then
        printf '%-28s OK       (%d s)\n' "$executable" "$elapsed"
    else
        printf '%-28s FAILED   (%d s, exit %d)\n' \
            "$executable" "$elapsed" "$status"

        failed_executables+=("$executable")
    fi
done

{
    echo "============================================================"
    echo "Sweep finished: $(date '+%Y-%m-%d %H:%M:%S %Z')"

    if ((${#failed_executables[@]} == 0)); then
        echo "Result: all ${#executables[@]} programs completed successfully."
    else
        echo "Result: ${#failed_executables[@]} program(s) failed:"
        printf '  %s\n' "${failed_executables[@]}"
    fi

    echo "============================================================"
} >> "$sweep_log"

echo

if ((${#failed_executables[@]} == 0)); then
    echo "All ${#executables[@]} tests completed successfully."
    echo "Detailed output: $sweep_log"
    exit 0
else
    echo "${#failed_executables[@]} of ${#executables[@]} tests failed." >&2
    echo "Detailed output: $sweep_log" >&2
    exit 1
fi
