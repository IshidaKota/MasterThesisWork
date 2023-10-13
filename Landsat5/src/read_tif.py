from skimage import io
import matplotlib.pyplot as plt
import numpy as np

def plot(img):
    fig,ax = plt.subplots()

    X,Y = np.meshgrid(img[:,0,0],img[0,:,0])

    ax.pcolormesh(X,Y,img[:,:,0].T)

    fig.savefig("../output/test_band1.png")
#path = "../19880330/LT05_L1TP_024033_19880330_20200917_02_T1_B1.TIF" # パスを格納
path = "../19880330/test.tif"
img = io.imread(path) 
print(img.shape) # 形の確認 #7021*7741 #numpy.ndarray #11135,5567,6
x,y,sr = img.shape
plot(img)




