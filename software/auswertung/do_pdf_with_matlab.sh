#!/bin/sh
echo "starting matlab and gerating figure. This may take some seconds..."
matlab -r simm -nodesktop -nosplash -r "pdffilename='ref'; create_pdf; exit" >> /dev/null
echo "...done!"
