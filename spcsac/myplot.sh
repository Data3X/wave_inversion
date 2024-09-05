#!/usr/bin/env -S bash -e
# GMT modern mode bash template
# Date:    2024-04-02T17:21:06
# User:    syr
# Purpose: Purpose of this script
# Set a unique session name
gmt begin map1 png
    gmt coast -R116.2/116.55/39.75/40.05 -JM12c -Baf -BWSne -W2p -A1000 -Glightbrown -Sazure1 --FORMAT_GEO_MAP=dddF
    gmt inset begin -DjBL+w3c/3.6c+o0.1c -F+gwhite+p1p
        gmt coast -R114/117/39/41 -JM? -EJP+glightbrown+p0.2p -A10000
        # 使用 -Sr+s 绘制矩形区域
        echo 116.2 39.75 116.55 40.05 | gmt plot -Sr+s -W1p,blue
    gmt inset end
gmt end show
