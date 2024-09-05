from obspy.imaging.beachball import beachball

# M11,M22,M33,M12,M13,M23
mt = [3.82, 3.39, -7.21, -4.55, -1.04, 5.04]
beachball(mt, size=200, linewidth=2, facecolor='b')

