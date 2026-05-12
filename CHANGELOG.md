# RELEASE NOTES

## Version 2.0.0 (Current)

### Features
- **Genome Alignment**: Replaced Bowtie2 transcriptome alignment with HISAT2 genome alignment.
- **BigWig Generation**:
  - Ribo-Seq: Generates strand-specific BigWig files (requires P-site correction).
  - RNA-Seq: Generates strand-specific BigWig files representing read coverage.
- **RNA-Seq Processing**: Full support for RNA-Seq data with genome alignment and deduplication options.

### Changes
- **Pipeline Logic**:
  - Removed `.ribo` file creation entirely. Genome-aligned reads cannot be represented in the transcriptome-coordinate `.ribo` format; this is now a pure genome-alignment pipeline producing bigWig and dedup BAM/BED artifacts.
  - Removed metadata-embedding flags (`do_metadata`, `input.metadata`, `root_meta`) — these only fed `.ribo` creation.
- **Environment**:
  - Renamed Conda environment to `ribo_genome`.
  - Docker image updated to `danielnguyener/riboflow`.

---

## Version 0.0.1
### Features
- Added UMI support
- Basic UMI extraction and deduplication

---

## Version 0.0.0

### Initial Release
- Initial RiboFlow pipeline release
- Transcriptome alignment support
- Basic deduplication methods
- RIBO file format support
