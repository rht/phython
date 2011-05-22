#http://zerokspot.com/weblog/2007/09/24/automating-stuff-with-scons/

import os
env=Environment(ENV=os.environ)
env.PDF(target = 'statmech/statmech.pdf', source = 'statmech/statmech.tex')
#Depends(pdfOutput,Split('presentation.tex bibliography.bib'))

