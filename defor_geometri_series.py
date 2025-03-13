import os
import numpy as np
from osgeo import gdal
from PyQt5.QtCore import QObject, pyqtSignal
import shutil
import matplotlib.pyplot as plt
# GDAL exceptions
gdal.UseExceptions()

#in_fn = fr"Y:\TGC-01530 Consejo COCOMACIA - TBD-CO-Choco\4_Risk\3_LULC_input\DF\T1T2_deforestation.tif"
#NRT = 3683
#n_classes = 30
def geometric_classification(self, in_fn, NRT, n_classes):
    '''
    geometric classification
    :param in_fn: map of distance from the forest eddge
    :param NRT:Negligible Risk Threshold
    :param n_classes:number of classes
    :return: mask_arr: result array with mask larger than NRT
    '''
    # Convert in_fn to NumPy array
    # Set up a GDAL dataset
    in_ds = gdal.Open(in_fn)
    # Set up a GDAL band
    in_band = in_ds.GetRasterBand(1)
    # Create Numpy Array
    arr = in_band.ReadAsArray()

    # The lower limit of the highest class = spatial resolution (the minimum distance possible without being in non-forest)
    LL = int(in_ds.GetGeoTransform()[1])

    self.progress_updated.emit(10)
    # The upper limit of the lowest class = the Negligible Risk Threshold
    UL = NRT = int(NRT)
    n_classes = int(n_classes)

    # Calculate common ratio(r)=(LLmax/LLmin)^1/n_classes
    r = np.power(LL / UL, 1/n_classes)

    # Create 2D class_array for the areas within the NRT
    class_array = np.array([[i, i + 1] for i in range(n_classes)])

    # Calculate UL and LL value for the areas within the NRT
    x= np.power(r, class_array)
    risk_class=np.multiply(UL,x)

    self.progress_updated.emit(20)
    # Create mask: areas beyond the NRT, assign class 1
    mask_arr=arr
    mask_arr[arr >= NRT] = 1

    self.progress_updated.emit(30)
    # Use boolean indexing to reclassification mask_arr value >= LL into risk_class
    # (e.g., if n_class is 29, class the areas within the NRT into class 2 to 30)
    # Set the progress_updated.emit() outside the loop to fasten the process
    mask_arr[(risk_class[0][0] > mask_arr) & (mask_arr >= risk_class[0][1])] = 2
    
    mask_arr[(risk_class[1][0] > mask_arr) & (mask_arr >= risk_class[1][1])] = 3
    mask_arr[(risk_class[2][0] > mask_arr) & (mask_arr >= risk_class[2][1])] = 4
    mask_arr[(risk_class[3][0] > mask_arr) & (mask_arr >= risk_class[3][1])] = 5
    mask_arr[(risk_class[4][0] > mask_arr) & (mask_arr >= risk_class[4][1])] = 6
    self.progress_updated.emit(40)
    mask_arr[(risk_class[5][0] > mask_arr) & (mask_arr >= risk_class[5][1])] = 7
    mask_arr[(risk_class[6][0] > mask_arr) & (mask_arr >= risk_class[6][1])] = 8
    mask_arr[(risk_class[7][0] > mask_arr) & (mask_arr >= risk_class[7][1])] = 9
    mask_arr[(risk_class[8][0] > mask_arr) & (mask_arr >= risk_class[8][1])] = 10
    mask_arr[(risk_class[9][0] > mask_arr) & (mask_arr >= risk_class[9][1])] = 11
    self.progress_updated.emit(50)
    mask_arr[(risk_class[10][0] > mask_arr) & (mask_arr >= risk_class[10][1])] = 12
    mask_arr[(risk_class[11][0] > mask_arr) & (mask_arr >= risk_class[11][1])] = 13
    mask_arr[(risk_class[12][0] > mask_arr) & (mask_arr >= risk_class[12][1])] = 14
    mask_arr[(risk_class[13][0] > mask_arr) & (mask_arr >= risk_class[13][1])] = 15
    mask_arr[(risk_class[14][0] > mask_arr) & (mask_arr >= risk_class[14][1])] = 16
    self.progress_updated.emit(60)
    mask_arr[(risk_class[15][0] > mask_arr) & (mask_arr >= risk_class[15][1])] = 17
    mask_arr[(risk_class[16][0] > mask_arr) & (mask_arr >= risk_class[16][1])] = 18
    mask_arr[(risk_class[17][0] > mask_arr) & (mask_arr >= risk_class[17][1])] = 19
    mask_arr[(risk_class[18][0] > mask_arr) & (mask_arr >= risk_class[18][1])] = 20
    mask_arr[(risk_class[19][0] > mask_arr) & (mask_arr >= risk_class[19][1])] = 21
    self.progress_updated.emit(70)
    mask_arr[(risk_class[20][0] > mask_arr) & (mask_arr >= risk_class[20][1])] = 22
    mask_arr[(risk_class[21][0] > mask_arr) & (mask_arr >= risk_class[21][1])] = 23
    mask_arr[(risk_class[22][0] > mask_arr) & (mask_arr >= risk_class[22][1])] = 24
    mask_arr[(risk_class[23][0] > mask_arr) & (mask_arr >= risk_class[23][1])] = 25
    mask_arr[(risk_class[24][0] > mask_arr) & (mask_arr >= risk_class[24][1])] = 26
    self.progress_updated.emit(80)
    mask_arr[(risk_class[25][0] > mask_arr) & (mask_arr >= risk_class[25][1])] = 27
    mask_arr[(risk_class[26][0] > mask_arr) & (mask_arr >= risk_class[26][1])] = 28
    mask_arr[(risk_class[27][0] > mask_arr) & (mask_arr >= risk_class[27][1])] = 29
    mask_arr[(risk_class[28][0] > mask_arr) & (mask_arr >= risk_class[28][1])] = 30
    self.progress_updated.emit(90)
    return mask_arr