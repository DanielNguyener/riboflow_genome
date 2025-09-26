# Modified RiboFlow Commands (rfc)

This is a modified version of the `rfcommands` package with additional features for genome alignment support in RiboFlow pipelines.

## New Features

### Added Commands
- **`rfc hisat2-log-to-csv`** - Convert HISAT2 alignment logs to CSV format
- **`rfc merge-hisat2-logs`** - Merge multiple HISAT2 log files

### Enhanced Commands
- **`rfc compile-step-stats`** - Added `--label-prefix` option for genome alignment stats
- **`rfc stats-percentage`** - Added `--label-prefix` option for genome alignment stats

## Installation

```bash
pip install -e .
```

## Verification

```bash
# Check that new commands are available
rfc --help | grep -E "(hisat2|merge-hisat2)"

# Test enhanced functionality
rfc compile-step-stats --help | grep "label-prefix"
rfc stats-percentage --help | grep "label-prefix"
```

## Compatibility

- **Python**: 3.7+
- **Dependencies**: click, pandas, numpy
- **Replaces**: Standard rfcommands package

## Usage Examples

```bash
# Convert HISAT2 log to CSV
rfc hisat2-log-to-csv -l alignment.log -n sample1 -p genome -o stats.csv

# Merge multiple HISAT2 logs
rfc merge-hisat2-logs -o merged.log log1.txt log2.txt log3.txt

# Compile stats with genome label prefix
rfc compile-step-stats --label-prefix genome input_dir output.csv
```