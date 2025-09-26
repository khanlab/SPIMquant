# Troubleshooting

This guide covers common issues encountered when using SPIMquant and provides solutions for typical problems.

## Installation Issues

### Missing Dependencies

**Problem:** `ModuleNotFoundError: No module named 'snakebids'`

**Solutions:**
```bash
# Install missing dependencies
pip install snakebids

# Or reinstall SPIMquant
pip install --upgrade git+https://github.com/khanlab/spimquant.git
```

**Problem:** `CommandNotFoundError: c3d not found`

**Solutions:**
```bash
# Use containers (recommended)
spimquant /data /output participant --use-apptainer

# Or install ITK-SNAP (includes c3d)
# On Ubuntu/Debian:
sudo apt install itksnap

# On macOS:
brew install --cask itk-snap
```

### Container Issues

**Problem:** `Container not found` or `Singularity/Apptainer not available`

**Solutions:**
```bash
# Install Apptainer (recommended)
# On Ubuntu 22.04+:
sudo apt install apptainer

# Alternative: use conda environments
spimquant /data /output participant --use-conda

# Or install dependencies natively
pip install greedyreg antspyx
```

## Data Input Issues

### BIDS Validation Errors

**Problem:** `BIDS validation failed`

**Solutions:**
```bash
# Skip validation temporarily
spimquant /data /output participant --skip-bids-validation

# Fix BIDS structure
# Ensure proper naming: sub-001_sample-brain_stain-PI_SPIM.ome.zarr
# Include dataset_description.json
```

**Example BIDS structure:**
```
dataset/
├── dataset_description.json
├── participants.tsv
└── sub-001/
    └── micr/
        ├── sub-001_sample-brain_stain-PI_SPIM.ome.zarr/
        └── sub-001_sample-brain_stain-abeta_SPIM.ome.zarr/
```

### OME-Zarr Format Issues

**Problem:** `Unable to read OME-Zarr file`

**Solutions:**
```python
# Check file structure
import zarr
store = zarr.open('file.ome.zarr', mode='r')
print(store.tree())

# Verify OME-Zarr compliance
from zarrnii import ZarrNii
znii = ZarrNii.from_ome_zarr('file.ome.zarr')
print(f"Channels: {znii.list_channels()}")
```

**Problem:** `Extension not recognized`

**Solutions:**
```bash
# For zip archives
spimquant /data /output participant --filter-spim extension='ome.zarr.zip'

# For directories  
spimquant /data /output participant --filter-spim extension='ome.zarr'
```

## Memory and Performance Issues

### Out of Memory Errors

**Problem:** `MemoryError` or `Killed` during processing

**Solutions:**
```bash
# Reduce resolution levels
spimquant /data /output participant --registration_level 6 --segmentation_level 3

# Limit parallel jobs
spimquant /data /output participant --cores 4 --resources mem_mb=8000

# Use more aggressive downsampling
spimquant /data /output participant --sloppy
```

**Memory usage guidelines:**
- Registration level 3: ~16GB per job
- Registration level 5: ~4GB per job  
- Registration level 6: ~2GB per job

### Slow Processing

**Problem:** Processing takes too long

**Solutions:**
```bash
# Use lower resolution for testing
spimquant /data /output participant --sloppy --registration_level 6

# Increase parallelization
spimquant /data /output participant --cores all

# Process subset of subjects
spimquant /data /output participant --participant_label 001 002
```

**Performance optimization:**
- Use SSD storage for work directory
- Ensure adequate RAM (2-4GB per core)
- Consider cloud instances with high memory

## Registration Issues

### Poor Registration Quality

**Problem:** Images poorly aligned to template

**Diagnostic steps:**
```bash
# Visual inspection
itksnap output/tpl-ABAv3/tpl-ABAv3_anat.nii.gz \
        output/sub-001/micr/sub-001_space-ABAv3_stain-PI_SPIM.nii

# Check registration stain
spimquant /data /output participant --stains_for_reg PI YOPRO AutoF
```

**Common causes and solutions:**

1. **Wrong stain for registration:**
   ```bash
   # Use nuclear/structural stains
   --stains_for_reg DAPI PI Hoechst AutoF
   ```

2. **Template mismatch:**
   ```bash
   # Try different template
   --template gubra  # Often better for lightsheet
   ```

3. **Intensity issues:**
   ```bash
   # Improve bias correction
   --correction_method n4
   ```

### Registration Failures

**Problem:** `Registration failed` error messages

**Solutions:**
```bash
# Check input data quality
# Ensure sufficient contrast
# Verify image orientation

# Try more robust registration
--registration_level 5  # Start with more downsampling

# Use hemisphere cropping if appropriate
--template_crop left
```

## Segmentation Issues

### Over-segmentation

**Problem:** Too much signal detected (high field fractions)

**Solutions:**
```bash
# Increase threshold
--seg_method threshold --seg_threshold 90

# Use more conservative Otsu
--seg_method otsu+k4i3  # More classes, higher class
```

### Under-segmentation

**Problem:** Too little signal detected (low field fractions)

**Solutions:**
```bash
# Decrease threshold
--seg_method threshold --seg_threshold 60

# Use more sensitive Otsu
--seg_method otsu+k3i1  # Lower class selection
```

### Segmentation Quality Check

```python
import nibabel as nib
import matplotlib.pyplot as plt

# Load original and segmented data
original = nib.load('sub-001_space-ABAv3_stain-abeta_SPIM.nii').get_fdata()
fieldfrac = nib.load('sub-001_seg-roi82_stain-abeta_fieldfrac.nii').get_fdata()

# Check segmentation quality
slice_idx = original.shape[2] // 2
plt.figure(figsize=(12, 4))

plt.subplot(1, 3, 1)
plt.imshow(original[:, :, slice_idx], cmap='gray')
plt.title('Original')

plt.subplot(1, 3, 2)
plt.imshow(fieldfrac[:, :, slice_idx], cmap='hot', vmin=0, vmax=1)
plt.title('Field Fraction')

plt.subplot(1, 3, 3)
plt.imshow(original[:, :, slice_idx], cmap='gray', alpha=0.7)
plt.imshow(fieldfrac[:, :, slice_idx], cmap='hot', alpha=0.3, vmin=0, vmax=1)
plt.title('Overlay')

plt.tight_layout()
plt.show()
```

## Cloud Storage Issues

### Authentication Problems

**Problem:** `Access denied` for cloud storage

**Amazon S3:**
```bash
# Set credentials
export AWS_ACCESS_KEY_ID=your_key
export AWS_SECRET_ACCESS_KEY=your_secret
export AWS_DEFAULT_REGION=us-east-1

# Or use profile
aws configure
export AWS_PROFILE=your_profile
```

**Google Cloud:**
```bash
# Service account
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json

# Or user credentials
gcloud auth application-default login
```

### Slow Cloud Access

**Problem:** Slow data transfer from cloud

**Solutions:**
```bash
# Increase parallelization
--jobs 20

# Use regional compute near data
# Ensure network bandwidth

# Cache frequently accessed data locally
--cache-storage /local/cache
```

## Output Issues

### Missing Output Files

**Problem:** Expected output files not generated

**Diagnostic steps:**
```bash
# Check for errors in logs
find output/logs -name "*.log" -exec grep -l "ERROR" {} \;

# Verify workflow completion
snakemake --summary

# Check specific rule status
snakemake --detailed-summary
```

### Corrupted Output

**Problem:** Output files appear corrupted or empty

**Solutions:**
```bash
# Clean and restart
rm -rf work/
spimquant /data /output participant --rerun-incomplete

# Check disk space
df -h

# Verify file integrity
file output/sub-001/micr/*.nii
```

## Workflow Issues

### Snakemake Errors

**Problem:** `WorkflowError` or `RuleException`

**Solutions:**
```bash
# Get detailed error information
spimquant /data /output participant --verbose

# Clean temporary files
rm -rf work/
spimquant /data /output participant --rerun-incomplete

# Debug specific rule
spimquant /data /output participant --debug-dag
```

### Lock File Issues

**Problem:** `Directory cannot be locked`

**Solutions:**
```bash
# Remove stale locks
find output/ -name ".snakemake_timestamp" -delete
find work/ -name ".snakemake_timestamp" -delete

# Or force unlock
snakemake --unlock
```

## Performance Monitoring

### Resource Usage

**Monitor resource consumption:**
```bash
# During processing
htop
free -h
df -h

# Post-processing analysis
cat output/logs/*/benchmark.txt
```

**Optimize resource allocation:**
```bash
# Based on monitoring, adjust:
--cores 8 --resources mem_mb=16000

# For cluster environments:
--cluster-config cluster.yaml
```

## Getting Help

### Log Analysis

**Key log locations:**
```
output/spimquant/logs/
├── import/          # Data loading issues
├── templatereg/     # Registration problems  
├── segmentation/    # Segmentation errors
└── quantification/  # Analysis issues
```

**Useful log commands:**
```bash
# Find errors
grep -r "ERROR" output/logs/

# Check recent failures
find output/logs -name "*.log" -mtime -1 -exec grep -l "fail\|error\|ERROR" {} \;

# Monitor progress
tail -f output/logs/templatereg/*.log
```

### Report Generation

**Generate diagnostic report:**
```bash
spimquant /data /output participant --report diagnostic.html --verbose
```

### Community Support

**When seeking help, include:**

1. **Command used:**
   ```bash
   spimquant /data /output participant --template ABAv3 --cores 8
   ```

2. **Error message:**
   ```
   Full error traceback from terminal
   ```

3. **System information:**
   ```bash
   python --version
   pip list | grep spim
   free -h
   df -h
   ```

4. **Data description:**
   - Number of subjects
   - Stains available
   - File sizes
   - Data source/format

**Resources:**
- GitHub Issues: https://github.com/khanlab/SPIMquant/issues
- Documentation: https://spimquant.readthedocs.io
- BIDS Specification: https://bids-specification.readthedocs.io