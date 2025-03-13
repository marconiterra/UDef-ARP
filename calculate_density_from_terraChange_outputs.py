# edited allocation tool
import os
import numpy as np
import pandas as pd
from osgeo import gdal
from PyQt5.QtCore import QObject, pyqtSignal
import shutil
import matplotlib.pyplot as plt
def set_working_directory(directory):
        '''
        Set up the working directory
        :param directory: your local directory with all dat files
        '''
        data_folder = directory
        os.chdir(data_folder)

###Step1 Create the Fitting Modeling Region Map###
def image_to_array(image):
        # Set up a GDAL dataset
        in_ds = gdal.Open(image)
        # Set up a GDAL band
        in_band = in_ds.GetRasterBand(1)
        # Create Numpy Array1
        arr = in_band.ReadAsArray()
        return arr

def array_to_image(in_fn, out_fn, data, data_type, nodata=None):
        '''
                Create image from array
                :param in_fn: datasource to copy projection and geotransform from
                :param out_fn:path to the file to create
                :param data:the NumPy array
                :param data_type:output data type
                :param nodata:optional NoData value
                :return:
        '''
        in_ds = gdal.Open(in_fn)
        output_format = out_fn.split('.')[-1].upper()
        if (output_format == 'TIF'):
                output_format = 'GTIFF'
        elif (output_format == 'RST'):
                output_format = 'rst'
        driver = gdal.GetDriverByName(output_format)
        out_ds = driver.Create(out_fn, in_ds.RasterXSize, in_ds.RasterYSize, 1, data_type, options=["BigTIFF=YES"])
        out_band = out_ds.GetRasterBand(1)
        out_ds.SetGeoTransform(in_ds.GetGeoTransform())
        out_ds.SetProjection(in_ds.GetProjection().encode('utf-8', 'backslashreplace').decode('utf-8'))
        if nodata is not None:
                out_band.SetNoDataValue(nodata)
        out_band.WriteArray(data)
        out_band.FlushCache()
        out_ds.FlushCache()
        return

def tabulation_bin_id_HRP(risk30_hrp, municipality, out_fn1):
        """
        This function is to create fitting modeling region array(tabulation_bin_id_masked)
        and fitting modeling region map(tabulation_bin_image)
        :param risk30_hrp: The 30-class vulnerability map for the CAL/HRP
        :param municipality: Subdivision image
        :param out_fn1: user input
        :return: tabulation_bin_id_masked: tabulation bin id array in CAL/HRP
        """
        # Convert risk30_hrp to NumPy array
        arr1 = image_to_array(risk30_hrp)

        # Convert municipality to NumPy array2
        arr2 = image_to_array(municipality)

        # Create a mask where the risk30_hrp value larger than 1 reclassed into 1
        mask_arr_HRP = np.where(arr1 > 0,1, arr1)

        # Calculate tabulation bin id with mask
        tabulation_bin_id_masked = np.add(arr1*1000, arr2) * mask_arr_HRP

        # Convert the array to signed 16-bit integer (int16) data type
        tabulation_bin_id_masked = tabulation_bin_id_masked.astype(np.int16)

        # Create the final image using tabulation_bin_image function
        array_to_image(risk30_hrp, out_fn1, tabulation_bin_id_masked,
                                     gdal.GDT_Int16, -1)
        return tabulation_bin_id_masked

def tabulation_bin_id_VP (risk30_vp, municipality, out_fn1):
        """
        This function is to create modeling region array(tabulation_bin_id_VP_masked)
        and modeling region map(tabulation_bin_image_vp)
        :param risk30_vp: The 30-class vulnerability map for the CNF/VP
        :param municipality: Subdivision image
        :param out_fn1: user input
        :return: tabulation_bin_id_VP_masked: tabulation bin id array in CNF/VP
        """

        # Convert municipality and risk30_vp to NumPy array
        arr2 = image_to_array(municipality)
        arr4 = image_to_array(risk30_vp)

        # Create a mask where the risk30_hrp value larger than 1 reclassed into 1
        mask_arr_VP = np.where(arr4 > 0, 1, arr4)

        # Calculate tabulation bin id of CNF/VP with mask
        tabulation_bin_id_VP_masked = np.add(arr4 * 1000, arr2) * mask_arr_VP

        # Convert the array to signed 16-bit integer (int16) data type
        tabulation_bin_id_VP_masked = tabulation_bin_id_VP_masked.astype(np.int16)
        #Write image to disk
        array_to_image(risk30_vp, out_fn1, tabulation_bin_id_VP_masked,
                                     gdal.GDT_Int16, -1)

        return tabulation_bin_id_VP_masked

def calculate_prediction_density_arr(risk30_vp, tabulation_bin_id_VP_masked, csv):
        '''
        Calculate the prediction density arry
        :param tabulation_bin_id_VP_masked: array for tabulation bin id in CNF/VP
        :param csv: relative frequency table
        :param risk30_vp: the 30-class vulnerability map for the CNF/VP
        :return: prediction_density_arr: modeled deforestation (MD)
        '''
        # Read Relative Frequency csv file
        merged_df=pd.read_csv(csv)

        # Insert index=0 row into first row of merged_df DataFrame
        new_row = pd.DataFrame({'ID': [0], 'Total Deforestation(pixel)': [0], 'Area of the Bin(pixel)': [0],
                                'Average Deforestation(pixel)': [0]})
        merged_df = pd.concat([new_row, merged_df]).reset_index(drop=True)

        # Using numpy.searchsorted() to assign values to 'id'
        df_sorted = merged_df.sort_values('ID')
        sorted_indices = df_sorted['ID'].searchsorted(tabulation_bin_id_VP_masked)
        relative_frequency_arr = tabulation_bin_id_VP_masked[:] = df_sorted['Average Deforestation(pixel)'].values[sorted_indices]

        # Calculate areal_resolution_of_map_pixels
        in_ds4 = gdal.Open(risk30_vp)
        P1 = in_ds4.GetGeoTransform()[1]
        P2 = abs(in_ds4.GetGeoTransform()[5])
        areal_resolution_of_map_pixels = P1 * P2 / 10000

        # Relative_frequency multiplied by the areal resolution of the map pixels to express the probabilities as densities
        prediction_density_arr=relative_frequency_arr * areal_resolution_of_map_pixels

        return prediction_density_arr

def replace_ref_system(in_fn, out_fn):
        '''
         RST raster format: correct reference system name in rdc file
         :param in_fn: datasource to copy correct projection name
         :param out_fn: rst raster file
        '''
        if out_fn.split('.')[-1] == 'rst':
            read_file_name, _ = os.path.splitext(in_fn)
            write_file_name, _ = os.path.splitext(out_fn)
            temp_file_path = 'rdc_temp.rdc'

            with open(read_file_name + '.rdc', 'r') as read_file:
                for line in read_file:
                    if line.startswith("ref. system :"):
                        correct_name=line
                        break

            if correct_name:
                with open(write_file_name + '.rdc', 'r') as read_file, open(temp_file_path, 'w') as write_file:
                    for line in read_file:
                        if line.startswith("ref. system :"):
                            write_file.write(correct_name)
                        else:
                            write_file.write(line)

                # Move the temp file to replace the original
                shutil.move(temp_file_path, write_file_name + '.rdc')

def calculate_adjustment_ratio(prediction_density_arr, expected_deforestation):
        '''
        Calculate the Adjustment Ratio (AR) in VP
        :param prediction_density_arr: modeled deforestation (MD)
        :param expected_deforestation: user input
        :return: AR
        '''

        # Sum up the pixels in the prediction density map. This is the modeled deforestation (MD).
        MD = np.sum(prediction_density_arr)

        # AR = ED / MD
        AR = expected_deforestation / MD
        return AR



def check_modeling_region_ids(csv, out_fn):
        '''
        Check modeling region IDs present in the prediction stage but absent in the fitting stage.
        :param csv: csv file of relative frequency in fitiing stage
        :param out_fn: modeling region image in prediction stage
        :return: id_difference: A set of modeling region IDs np array that exist only in the prediction stage
        '''
        fit_model_region_id = pd.read_csv(csv)['ID'].to_numpy()
        pre_model_region_arr = image_to_array(out_fn)
        pre_model_region_id = np.unique(pre_model_region_arr[pre_model_region_arr != 0])
        id_difference = np.setdiff1d(pre_model_region_id, fit_model_region_id)

        return id_difference


def adjusted_prediction_density_map_annual (prediction_density_arr, risk30_vp, AR, out_fn2, time):
        '''
        Create adjusted prediction density map for annual
        :param prediction_density_arr:modeled deforestation (MD)
        :param risk30_vp: risk30_vp image
        :param AR:Adjustment Ratio
        :param out_fn2: user input
        :return:
        '''

        # Calculate the maximum density
        # Calculate areal_resolution_of_map_pixels
        in_ds4 = gdal.Open(risk30_vp)
        P1 = in_ds4.GetGeoTransform()[1]
        P2 = abs(in_ds4.GetGeoTransform()[5])
        maximum_density = P1 * P2 / 10000

        # Adjusted_Prediction_Density_Map = AR x Prediction_Density _Map
        adjusted_prediction_density_arr=AR*prediction_density_arr

        # Reclassify all pixels greater than the maximum (e.g., 0.09) to be the maximum
        adjusted_prediction_density_arr[adjusted_prediction_density_arr > maximum_density] = maximum_density

        # Convert the result back to an annual rate by dividing by the number of years in the VP
        adjusted_prediction_density_arr_annual=adjusted_prediction_density_arr/time

        # Create imagery
        self.array_to_image(risk30_vp, out_fn2, adjusted_prediction_density_arr_annual, gdal.GDT_Float32, -1)

        return

def calculate_missing_bins_rf (id_difference, csv):
        '''
        If one or more empty bins are found, compute the jurisdiction-wide weighted average of relative frequencies for
        missing bins and create a new csv file
        :param csv: csv file of relative frequency in the fitting stage
        :param id_difference: A set of modeling region IDs np array that exist only in the prediction stage
        :return
        '''
        # Convert modeling region ids to vulnerability zone id
        df=pd.read_csv(csv)
        df['v_zone'] = (df['ID'] // 1000).astype(int)

        # Convert missing bin ids to vulnerability zone id
        missing_v_zone = [x // 1000 for x in id_difference]

        # Select rows
        filtered_df = df[df['v_zone'].isin(missing_v_zone)].copy()

        # Created new column
        filtered_df['Total Deforestation(pixel)'] = filtered_df['Area of the Bin(pixel)'] * filtered_df['Average Deforestation(pixel)']

        # Group by and sum area and weighted relative frequency
        aggregated_df = filtered_df.groupby('v_zone')[['Total Deforestation(pixel)', 'Area of the Bin(pixel)']].sum().reset_index()

        # Calculate Average Deforestation
        aggregated_df['Average Deforestation(pixel)']=aggregated_df['Total Deforestation(pixel)']/aggregated_df['Area of the Bin(pixel)']

        # Create a new dataframe id_difference_df
        id_difference_df = pd.DataFrame(id_difference, columns=['ID'])
        id_difference_df['v_zone'] = missing_v_zone

        # Create missing bins dataframe by outer join aggregated_df and id_difference_df
        missing_bins_df = pd.merge(id_difference_df,aggregated_df , on='v_zone', how='outer')

        # Insert missing bins dataframe back to csv file
        df_new = pd.concat([df, missing_bins_df], ignore_index=True)

        # Sorting by column "ID"
        df_new=df_new.sort_values(by=['ID'], ascending=True)

        # Drop column 'v_zone'
        df_new=df_new.drop(['v_zone'], axis=1)

        # Copy the original csv file copy and rename it to csv_orig
        shutil.copyfile(csv, csv.split('.')[0] + '_orig' + '.csv')

        # Save the new result to csv
        df_new.to_csv(csv, index=False)


def adjusted_prediction_density_array (prediction_density_arr, risk30_vp, AR):
        '''
        Create adjusted prediction density array
        :param prediction_density_arr:modeled deforestation (MD)
        :param risk30_vp: risk30_vp image
        :param AR:Adjustment Ratio
        :return: adjusted_prediction_density_np_arr
        '''

        # Calculate the maximum density
        # Calculate areal_resolution_of_map_pixels
        in_ds4 = gdal.Open(risk30_vp)
        P1 = in_ds4.GetGeoTransform()[1]
        P2 = abs(in_ds4.GetGeoTransform()[5])
        maximum_density = P1 * P2 / 10000

        # Adjusted_Prediction_Density_Map = AR x Prediction_Density _Map
        adjusted_prediction_density_arr=AR*prediction_density_arr

        # Reclassify all pixels greater than the maximum (e.g., 0.09) to be the maximum
        adjusted_prediction_density_arr[adjusted_prediction_density_arr > maximum_density] = maximum_density

        return adjusted_prediction_density_arr


# open tabulation_bin_id as array
def tabulation_bin_id_HRP(risk30_hrp, municipality, out_fn1):
        """
        This function is to create fitting modeling region array(tabulation_bin_id_masked)
        and fitting modeling region map(tabulation_bin_image)
        :param risk30_hrp: The 30-class vulnerability map for the CAL/HRP
        :param municipality: Subdivision image
        :param out_fn1: user input
        :return: tabulation_bin_id_masked: tabulation bin id array in CAL/HRP
        """
        # Convert risk30_hrp to NumPy array
        arr1 = image_to_array(risk30_hrp)

        # Convert municipality to NumPy array2
        arr2 = image_to_array(municipality)

        # Create a mask where the risk30_hrp value larger than 1 reclassed into 1
        mask_arr_HRP = np.where(arr1 > 0,1, arr1)

        # Calculate tabulation bin id with mask
        tabulation_bin_id_masked = np.add(arr1*1000, arr2) * mask_arr_HRP

        # Convert the array to signed 16-bit integer (int16) data type
        tabulation_bin_id_masked = tabulation_bin_id_masked.astype(np.int16)

        # Create the final image using tabulation_bin_image function
        array_to_image(risk30_hrp, out_fn1, tabulation_bin_id_masked,
                                     gdal.GDT_Int16, -1)
        return tabulation_bin_id_masked

def create_relative_frequency_table(tabulation_bin_id_masked, deforestation_hrp, csv_name):
        """
        Create dataframe
        :param tabulation_bin_id_masked: array with id and total deforestation
        :param deforestation_hrp: Deforestation Map during the CAL/HRP
        :return: merged_df: relative frequency dataframe
        """
        # Calculate array area of the bin [integer] (in pixels) for Col3 using np.unique and counts function, excluding 0
        unique, counts = np.unique(tabulation_bin_id_masked[tabulation_bin_id_masked != 0], return_counts=True)
        # Convert to array
        arr_counts = np.asarray((unique, counts)).T

        # Calculate total deforestation within the bin [integer] array for Col2
        arr3 = image_to_array(deforestation_hrp)

        # deforestation_within_bin will have tabulation_bin_id value in deforestation pixel
        deforestation_within_bin = tabulation_bin_id_masked * arr3
        # Use np.unique to counts total deforestation in each bin
        unique1, counts1 = np.unique(deforestation_within_bin[deforestation_within_bin != 0], return_counts=True)


        # Convert to array
        arr_counts_deforestion = np.asarray((unique1, counts1)).T

        # Create pandas DataFrames
        df1 = pd.DataFrame(arr_counts_deforestion, columns=['ID', 'Total Deforestation(pixel)'])
        df2 = pd.DataFrame(arr_counts, columns=['ID', 'Area of the Bin(pixel)'])

        # Merge the two DataFrames based on the 'id' column using an outer join to include all rows from both DataFrames
        merged_df = pd.merge(df1, df2, on='ID', how='outer').fillna(0)

        # Calculate Average Deforestation by performing the division operation of col2 and col3 and add a new column to merged_df
        merged_df['Average Deforestation(pixel)'] = merged_df.iloc[:, 1].astype(float) / merged_df.iloc[:, 2].astype(float)

        # Sort the DataFrame based on the 'ID'
        merged_df = merged_df.sort_values(by='ID')

        # Reset the index to have consecutive integer indices
        merged_df = merged_df.reset_index(drop=True)

        csv_file_path = csv_name
        merged_df.to_csv(csv_file_path, index=False)

        return merged_df

def execute_workflow_vp(directory,max_iterations, csv, municipality, expected_deforestation, risk30_vp, out_fn1, out_fn2, time):
        '''
        Create workflow function for VP
        :param max_iterations: maximum number of iterations
        '''
        data_folder = set_working_directory(directory)
        tabulation_bin_id_VP_masked = tabulation_bin_id_VP(risk30_vp, municipality, out_fn1)
        replace_ref_system(municipality, out_fn1)

        # Check modeling region IDs present in the prediction stage but absent in the fitting stage
        id_difference = check_modeling_region_ids(csv, out_fn1)

        # If there are missing bins, calculate the relative frequency and create a new csv file
        if id_difference.size > 0:
                calculate_missing_bins_rf(id_difference, csv)


        prediction_density_arr = calculate_prediction_density_arr(risk30_vp, tabulation_bin_id_VP_masked,csv)
        AR = calculate_adjustment_ratio(prediction_density_arr, expected_deforestation)
        # Set a maximum number of iterations to avoid infinite loop
        iteration_count = 0
        new_prediction_density_arr = None

        # When AR > 1.00001 and iteration_count <= max_iterations, treat the result as new prediction density map to iterate the AR util AR is <=1.00001 or iteration_count = max_iterations
        while AR > 1.00001 and iteration_count <= max_iterations:
                new_prediction_density_arr = adjusted_prediction_density_array(prediction_density_arr, risk30_vp, AR)
                AR = calculate_adjustment_ratio(new_prediction_density_arr, expected_deforestation)
                iteration_count += 1
                # Emitting progress based on the current iteration_count and max_iterations
        if iteration_count <= int(max_iterations):
                selected_density_arr = new_prediction_density_arr if new_prediction_density_arr is not None else prediction_density_arr
                adjusted_prediction_density_map_annual(selected_density_arr, risk30_vp, AR, out_fn2, time)
                replace_ref_system(municipality, out_fn2)
        else:
                print("Maximum number of iterations reached. Please reset the maximum number of iterations.")


        return id_difference

def reclassify_raster(raster, mapping):
    reclassified = np.full(raster.shape, 'NF', dtype=object)  # Default to 'NF'
    for class_num, category in mapping.items():
        reclassified[raster == class_num] = category
    return reclassified

def get_deforestation_from_TerraChange(initial_year, final_year, in_path, out_path, class_definition_file):
        '''
        Get deforestation from TerraChange outputs
        :param initial_year: initial year
        :param final_year: final year
        :param in_path: input path
        :param out_path: output path
        :param class_definition_file: class definition file
        :return: deforestation_hrp
        '''
        list_rasters = [f for f in os.listdir(in_path) if f.endswith('.tif')]

        # Create a mapping from classNumber to categoryList ('F' or 'NF')
        class_definition = pd.read_csv(class_definition_file)
        
        class_mapping = dict(zip(class_definition['classNumber'], class_definition['categoryList']))
        list_rasters.sort()
        start_path = list_rasters[initial_year]
        end_path = list_rasters[final_year]

        # Open the start and end rasters
        start_ds =  image_to_array(fr"{in_path}//{start_path}")
        start_ds = reclassify_raster(start_ds, class_mapping)
        end_ds = image_to_array(fr"{in_path}//{end_path}")
        end_ds = reclassify_raster(end_ds, class_mapping)

        # Calculate the deforestation by subtracting the end raster from the start raster
        # Create a deforestation mask
        deforestation_mask = np.where((start_ds == 'F') & (end_ds == 'NF'), 1, 0)

        # Convert the mask to int8 data type to save memory
        deforestation_hrp = deforestation_mask.astype(np.int8)

        # Save the deforestation map
        array_to_image(fr"{in_path}//{start_path}", out_path, deforestation_hrp, gdal.GDT_Byte, 0)
        

        return deforestation_hrp

in_fn = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\2_LULC_input\DF\T1_distance.tif"
deforestation_hrp = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\2_LULC_input\DF\T1T2_deforestation.tif"
mask = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\2_jurisdiction_and_region\JNR_region_SF_IHR.tif"
out_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\3_output_CAL\v33"
def nrt_calculation(in_fn, deforestation_hrp, mask, out_path):
        '''
        NRT calculation
        :param in_fn: map of distance from the forest eddge in CAL
        :param deforestation_hrp:deforestation binary map in HRP
        :param mask: mask of the non-excluded jurisdiction (binary map)
        :return: NRT: Negligible Risk Threshold
        '''
        # Convert image to NumPy array
        distance_arr_cal = image_to_array(in_fn)
        deforestation_hrp_arr = image_to_array(deforestation_hrp)
        mask_arr = image_to_array(mask)

        # Mask the distance arr within deforstation pixel and study area
        distance_arr_masked=distance_arr_cal*mask_arr*deforestation_hrp_arr

        ## Calculate the histogram
        # Flatten the distance_arr_masked and expect 0 for np.histogram function
        # The np.histogram is computed over the flattened array
        distance_arr_masked_1d = distance_arr_masked.flatten()
        distance_arr_masked_1d = distance_arr_masked_1d[distance_arr_masked_1d != 0]

        ## Calculate the histogram
        # Set up bin width as spatial resolution
        in_ds = gdal.Open(in_fn)
        P = in_ds.GetGeoTransform()[1]
        bin_width =int(P)
        # Calculate the histogram
        hist, bin_edges = np.histogram(distance_arr_masked_1d, bins=np.arange(distance_arr_masked_1d.min(),
                                                                              distance_arr_masked_1d.max() + bin_width,
                                                                               bin_width))
        # plot the histogram
        plt.figure(figsize=(10, 6))
        plt.bar(bin_edges[:-1], hist, width=bin_width, align='edge')
        plt.xlabel('Distance from forest edge (m)')
        plt.ylabel('Frequency of deforestation (pixels)')
        plt.title('Histogram of Deforestation Distance from Forest Edge')
        # save plot as png
        plt.savefig(fr'{out_path}/histogram_deforestation_distance.png')
        plt.close()

        # Calculate the cumulative proportion
        # Normalize the histogram to get probability
        hist_normalized = hist / np.sum(hist)

        # Compute cumulative distribution
        cumulative_prop = np.cumsum(hist_normalized)

        # # Find the index cumulative proportion >= 0.995
        index_995 = np.argmax(cumulative_prop >= 0.995)
        
        # Get the bin edges for the NRT bin
        nrt_bin_start = bin_edges[index_995]
        nrt_bin_end = bin_edges[index_995 + 1]
        # Create a cumulative histogram
        plt.figure(figsize=(10, 6))
        
        # Plot bins before index_995 in blue
        plt.bar(bin_edges[:-1][:index_995], cumulative_prop[:index_995], 
                width=bin_width, align='edge', color='blue', alpha=0.7)
        
        # Plot bins after index_995 in green
        plt.bar(bin_edges[:-1][index_995:], cumulative_prop[index_995:], 
                width=bin_width, align='edge', color='green', alpha=0.7)
        
        # Plot vertical line at bin_edges[index_995]
        plt.axvline(x=bin_edges[index_995], color='red', linestyle='--', linewidth=2)
        
        plt.xlabel('Distance from forest edge (m)')
        plt.ylabel('Cumulative proportion')
        plt.title('Cumulative Histogram of Deforestation Distance from Forest Edge')
        plt.legend(['NRT threshold', 'Below NRT', 'Above NRT'])
        
        # Add text annotation for NRT value
        plt.text(bin_edges[index_995], 0.5, f'NRT: {bin_edges[index_995]:.2f}m', 
                 rotation=90, verticalalignment='center')
        
        # save plot as png
        plt.savefig(fr'{out_path}/cumulative_histogram_deforestation_distance.png')
        plt.close()

        # Calculate the average of the NRT bin
        NRT = int((nrt_bin_start + nrt_bin_end) / 2)
        return NRT



risk30_hrp = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\4_output_application_VB\ALT2_VP_Vulnerability.tif"
municipality = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\2_jurisdiction_and_region\JNR_Regions.tif"
class_definition_file = fr"Y:\\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\2_RS\4_Auxilary\class definition.csv"
initial_year = 0
final_year = 6
in_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\2_empirical_transition_potential\terraChange_outputs\VP\Predictions_Projected"
out_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\4_output_application_VP\VP_TC_deforestation.tif"
out_fn1 = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\4_output_application_VB\out_fn1.tif"
deforestation_hrp = get_deforestation_from_TerraChange(initial_year, final_year, in_path, out_path, class_definition_file)
csv_name = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\4_output_application_VB\data.csv"
#tabulation_bin_id_masked = tabulation_bin_id_HRP(risk30_hrp, municipality, out_fn1)
deforestation_vp = out_path
#create_relative_frequency_table(tabulation_bin_id_masked, deforestation_vp, csv_name)
import rasterio
def calculate_deforestation_density(deforestation_hrp, output_path):
    """
    Calculate deforestation density in 100m x 100m grid cells.
    
    Args:
    deforestation_hrp (str): Path to the input deforestation raster.
    output_path (str): Path to save the output density raster.
    
    Returns:
    None
    """
    # Open the deforestation raster
    with rasterio.open(deforestation_hrp) as src:
        deforestation = src.read(1)
        profile = src.profile.copy()
        
        # Calculate the number of pixels in a 100m x 100m grid
        pixels_per_grid = int(100 / src.res[0])
        
        # Pad the array if necessary to ensure it's divisible by pixels_per_grid
        pad_rows = (pixels_per_grid - deforestation.shape[0] % pixels_per_grid) % pixels_per_grid
        pad_cols = (pixels_per_grid - deforestation.shape[1] % pixels_per_grid) % pixels_per_grid
        deforestation_padded = np.pad(deforestation, ((0, pad_rows), (0, pad_cols)), mode='constant', constant_values=0)
        
        # Reshape the array to group pixels into 100m x 100m grids
        reshaped = deforestation_padded.reshape(
            deforestation_padded.shape[0] // pixels_per_grid, 
            pixels_per_grid, 
            deforestation_padded.shape[1] // pixels_per_grid, 
            pixels_per_grid
        )
        
        # Calculate the sum of deforestation pixels (value 1) in each grid
        deforestation_sum = np.sum(reshaped == 1, axis=(1, 3))
        
        # Calculate the total number of pixels in each grid
        total_pixels = pixels_per_grid * pixels_per_grid
        
        # Calculate density (sum of 1s over total pixels in the grid)
        density = deforestation_sum / total_pixels
        
        # Resize the density array to match the original raster dimensions
        density_resized = np.repeat(np.repeat(density, pixels_per_grid, axis=0), pixels_per_grid, axis=1)
        density_resized = density_resized[:deforestation.shape[0], :deforestation.shape[1]]
        
        # Update the profile for the output raster
        profile.update(dtype=rasterio.float32, count=1, compress='lzw')
        
        # Write the output raster
        with rasterio.open(output_density_path, 'w', **profile) as dst:
            dst.write(density_resized.astype(rasterio.float32), 1)

# Example usage
deforestation_hrp = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\4_output_application_HRP\HRP_TC_deforestation.tif"
output_density_path = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\4_output_application_HRP\HRP_TC_ETP.tif"
calculate_deforestation_density(deforestation_hrp, output_density_path)

