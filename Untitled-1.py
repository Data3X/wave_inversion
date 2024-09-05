from obspy import read
from matplotlib import pyplot as plt
singlechannel = read('D:\\research\\dsm\\DSMsynTI-mpi-master\\spcsac\\example\\XYZ.20090101.Rs')
data = singlechannel[0].data
length = len(data)
lwin = 400
swin = 10
sta = 0
lta = 0
STA_LTA = [None for i in range(lwin+swin)]
def func(start,t):
    return data[t+start]**2

for start in range(length-lwin-swin):
    for lcount in range(lwin):
        lta = lta + func(start,lcount)
    for scount in range(swin):
        sta = sta + func(start+lwin,scount)
    STA_LTA.append(sta*lwin/lta/swin)

plt.plot(STA_LTA)

singlechannel.plot()
#singlechannel[0].spectrogram()
singlechannel[0].spectrogram(log=True,samp_rate=12.8,per_lap=0.95,wlen=200,dbscale=True)
# for count in len(data):
#     data[count]
# 绘制单分量波形数据，默认大小为 800x250
# 添加 outfile 参数后图片保存到本地
#singlechannel.plot()