#gauge近傍の領域抽出
import matplotlib.pyplot as plt
def extract_gaugeRect(coordinate,bands2d,boxRange=0.001):
    bboxXMin,bboxXMax = coordinate[0]-boxRange,coordinate[0]+boxRange
    bboxYMin,bboxYMax = coordinate[1]-boxRange,coordinate[1]+boxRange
    print(f"xrange={bboxXMin,bboxXMax},yrange={bboxYMin,bboxYMax}")
    rect = bands2d.sel(x=slice(bboxXMin,bboxXMax),y=slice(bboxYMax,bboxYMin))
    rect.isel(band=0).plot()
    return rect


def extract_RiverReflectance(rect,title="",figname=""):
    bandlen,ylen,xlen =rect.shape
    #region = region.dropna(dim="band",how="all")
    #print(region)
    fig,ax = plt.subplots()
    wavelength = [485,565,660,830]
    spectrumList = []
    for i in range(1,ylen-1):
        for j in range(1,xlen-1):
            if rect[:,i-2:i+2,j-2:j+2].all() != 0:  
            #if region[0,i-1,j-1] != 0 and region[0,i+1,j-1] != 0 and region[0,i-1,j+1] != 0 and region[0,i-1,j+1] != 0:
                ax.plot(wavelength,rect[:4,i,j])
                spectrumList.append(list(rect[:4,i,j]))
    
    #modify
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Reflectance (sr^-1)")
    plt.show()

    return spectrumList

def aveReflectance(spectrumList):
    #average reflectance
    band1,band2,band3,band4 = 0,0,0,0
    for spectrum in spectrumList:
        band1 += spectrum[0]
        band2 += spectrum[1]
        band3 += spectrum[2]
        band4 += spectrum[3]
    band1 /= len(spectrumList)
    band2 /= len(spectrumList)
    band3 /= len(spectrumList)
    band4 /= len(spectrumList)

    #plot
    fig,ax = plt.subplots()
    wavelength = [485,565,660,830]
    ax.plot(wavelength,[band1,band2,band3,band4])
    plt.show()
    return [band1,band2,band3,band4]