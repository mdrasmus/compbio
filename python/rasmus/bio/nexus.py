

class Nexus (object):
    
    def __init__(self, filename=None, mod="w"):
        self.out = None
        
        if filename is not None and mod == "w":
            self.write_header(filename)
            
    
    def write_header(self, filename):
        if isinstance(filename, basestring):
            self.out = open(filename, "w")
        else:
            self.out = filename
        
        self.out.write("#NEXUS\n")
    
    
    def write_matrix(self, names, seqs, format, seqwidth=1000):
        
        if format == "pep":
            format = "protein"
        
        self.out.write("""begin data;
dimensions ntax=%d nchar=%d;
format datatype=%s interleave=yes gap=-;
matrix

""" % (len(seqs), len(seqs[0]), format))
    
        totalsize = len(seqs[0])
        size = 0

        # write sequences
        while size < totalsize:
            for name, seq in zip(names, seqs):
                self.out.write("%s  %s\n" % (name, seq[size:size+seqwidth]))
            size += seqwidth

        self.out.write(";\nend;\n")


    def write(self, data):
        self.out.write(data)

    
    def close(self):
        self.out.close()
