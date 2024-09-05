#!/bin/bash
gmt begin GlobalMap
    gmt coast -JM12c -R100/150/10/50 -Ba -W0.5p -A10000
    gmt plot -Gred -Sa0.3c -W0.5p -l"Event" << EOF
145.59  18.83
EOF
echo 111.804 36.634 118.804 42.634 | gmt plot -Sr+s -W1p,blue
    #gmt inset begin -DjBL+w3c/3.6c+o0.1c -F+gwhite+p1p
    #    gmt coast -R111.804/118.804/36.634/42.634 -JM? -EJP+glightbrown+p0.2p -A10000
    #    echo 100 10 150 50 | gmt plot -Sr+s -W1p,blue
    #gmt inset end
gmt end show