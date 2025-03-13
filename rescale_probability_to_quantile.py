import numpy as np
import rasterio


def rescale_raster_to_quantiles(input_raster_path, output_raster_path, n_bins=30, nodata_value=0):
    """
    Rescale positive raster values into quantile distribution bins.
    Parameters:
    input_raster_path (str): Path to input raster file
    output_raster_path (str): Path to save output raster file
    n_bins (int): Number of quantile bins (default=30)
    nodata_value (float, optional): Value to treat as no data
    Returns:
    None: Saves the rescaled raster to output_raster_path
    """
    with rasterio.open(input_raster_path) as src:
        # Read the raster data
        raster_data = src.read(1)
        profile = src.profile.copy()
        # Create mask for valid data
        if nodata_value is None:
            nodata_value = src.nodata
        if nodata_value is not None:
            valid_mask = (raster_data != nodata_value) & (raster_data > 0)
        else:
            valid_mask = raster_data > 0
        # Get valid data for percentile calculation
        valid_data = raster_data[valid_mask]
        if len(valid_data) == 0:
            raise ValueError("No valid data found in raster")
        # Calculate percentile boundaries
        percentiles = np.percentile(valid_data, np.linspace(0, 100, n_bins+1))
        percentiles = np.unique(percentiles)
        # Handle edge case where there aren't enough unique value
        # Create output array
        output_data = np.zeros_like(raster_data)
        # Assign bins only to valid data
        output_data[valid_mask] = np.clip(
            np.digitize(raster_data[valid_mask], percentiles) - 1,
            0, n_bins-1
        )
        # Set nodata values
        if nodata_value is not None:
            output_data[~valid_mask] = nodata_value
        # Update the profile for output
        profile.update(
            dtype=rasterio.uint8 if n_bins <= 255 else rasterio.uint16,
            nodata=nodata_value
        )
        # rescale to 0-1
        output_data = output_data / (n_bins)
        # 
        # fill data originally nodata_value with 0
        output_data[raster_data == 0] = 0
        # Write the output raster
        profile.update(dtype=rasterio.float32)
        # remove the nodata value
        profile.update(nodata=None)
        with rasterio.open(output_raster_path, 'w', **profile) as dst:
            dst.write(output_data.astype('float32'), 1)

def get_quantile_stats(input_raster_path, n_bins=3000, nodata_value=None):
    """
    Get statistics about the quantile distribution of the raster.
    Parameters:
    input_raster_path (str): Path to input raster file
    n_bins (int): Number of quantile bins (default=30)
    nodata_value (float, optional): Value to treat as no data
    Returns:
    dict: Statistics about the quantile distribution
    """
    with rasterio.open(input_raster_path) as src:
        raster_data = src.read(1)
        if nodata_value is None:
            nodata_value = src.nodata
        if nodata_value is not None:
            valid_mask = (raster_data != nodata_value) & (raster_data > 0)
        else:
            valid_mask = raster_data > 0
        valid_data = raster_data[valid_mask]
        if len(valid_data) == 0:
            return None
        percentiles = np.percentile(valid_data, np.linspace(0, 100, n_bins+1))
        return {
            'bin_boundaries': percentiles,
            'min_value': valid_data.min(),
            'max_value': valid_data.max(),
            'mean_value': valid_data.mean(),
            'median_value': np.median(valid_data),
            'n_valid_pixels': len(valid_data)
        }
    

version= 'v40'
transition = "DG"
input_raster_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\{version}\3_LULC_input\{transition}\HRP_ETP.tif"
output_raster_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\{version}\3_LULC_input\{transition}\HRP_qETP.tif"
rescale_raster_to_quantiles(input_raster_path, output_raster_path, n_bins=30, nodata_value=0)

#input_raster_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\{version}\3_LULC_input\{transition}\VP_ETP.tif"
#output_raster_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\{version}\3_LULC_input\v38\{transition}\VP_qETP.tif"
#rescale_raster_to_quantiles(input_raster_path, output_raster_path, n_bins=30, nodata_value=0)

input_raster_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\{version}\3_LULC_input\{transition}\PRD_ETP.tif"
output_raster_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\{version}\3_LULC_input\{transition}\PRD_qETP.tif"
rescale_raster_to_quantiles(input_raster_path, output_raster_path, n_bins=30, nodata_value=0)

input_raster_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\{version}\3_LULC_input\{transition}\CAL_ETP.tif"
output_raster_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\{version}\3_LULC_input\{transition}\CAL_qETP.tif"
rescale_raster_to_quantiles(input_raster_path, output_raster_path, n_bins=30, nodata_value=0)