import matplotlib
# matplotlib.use( "agg" )

# attempt to import the PyPDF2 module, if it is not present then 
# the annotation function just does nothing
no_pdf_module = False
try:
    import PyPDF2 as pdf
except ImportError:
    print ("No PyPDF2 installed, figure annotations disabled")
    no_pdf_module = True
    pass


def fig_annotate(filename, params):
    '''
    Given a figure and a dictionary of values (like similation parameter values) 
    associated with that figure, write in the "Keyword" metadata of the PDF the parameter list.
    
    When conducting numerical simulations involving large numbers of parameters, it 
    can be difficult and confusing to keep track of which parameters are associated with which
    files. This function provides an additional check by allowing paramter values to be written directly 
    into the PDF file
    
    William Gilpin, 2015
    '''

    if no_pdf_module:
        pass
    else:
        inpdf = pdf.PdfFileReader(filename)
        
        mypage = inpdf.getPage(0)
        outpdf = pdf.PdfFileWriter()
        outpdf.addPage(mypage)
        metadata_string = ''
        for item in params:
            metadata_string += str(item) + '='+str(params[item])
            metadata_string += '    '
        outpdf.addMetadata({'/Keywords': metadata_string})


        outputStream = open(filename, 'wb')
        outpdf.write(outputStream)
    outputStream.close()
    
def savefig2(filename, params):
    '''
    A wrapper for matplotlib's savefig that includes metadata
    
    William Gilpin, 2015
    '''
    savefig(filename)
    fig_annotate(filename, params)
    
def fig_get_annotations():
    '''
    Given a figure with annotations written into the "keyword" field 
    using the format of fig_annotate.py, read those parameter values 
    and return a parameter dictionary
    
    NOT YET IMPLEMENTED
    '''
    
    
    
    return None