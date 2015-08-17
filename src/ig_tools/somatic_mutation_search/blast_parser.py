import os
import sys

class BlastBlock:
    cluster_id = -1
    matches = list()

    def __init__(self):
        self.cluster_id = -1
        self.matches = list()

    def ParseClusterId(self, id_str):
        splits = id_str.split('___')
        if len(splits) != 4:
            print("ERROR: Query string " + id_str + " is not correct")
            sys.exit(1)
        self.cluster_id = splits[1]

    def Empty(self):
        return self.cluster_id == -1 and len(self.matches) == 0

    def Print(self):
        print("Id: " + str(self.cluster_id))
        print(self.matches)

class BlastSettings:
    query_str = "Query= "
    align_str = "Sequences producing significant alignments"
    num_best_matches = 3
    

class BlastOutput:
    id_map = dict()
    blocks = list()

    def Add(self, blast_block):
        self.blocks.append(blast_block)
        self.id_map[blast_block.cluster_id] = len(self.blocks) - 1

    def size(self):
        return len(self.blocks)

    def Parse(self, output_fname):
        if not os.path.exists(output_fname):
            print("ERROR: File with BLAST output " + output_fname + " was not found")
            sys.exit(1)

        lines = open(output_fname, "r").readlines()
        blast_block = BlastBlock()
        start_align = False
        num_processed_matches = 0
        print(str(len(lines)) + " lines were extracted from " + output_fname)
        for l in lines:
            l = l.strip()
            if l[:len(BlastSettings.query_str)] == BlastSettings.query_str:
                if not blast_block.Empty():
                    self.Add(blast_block)
                    blast_block = BlastBlock()
                    start_align = False
                    num_processed_matches = 0
                blast_block.ParseClusterId(l)
            elif l[:len(BlastSettings.align_str)] == BlastSettings.align_str:
                start_align = True
            elif start_align and l != "" and num_processed_matches <= BlastSettings.num_best_matches:
                splits = l.split(' ')
                blast_block.matches.append(splits[0])
                num_processed_matches += 1
        self.Add(blast_block)    
        print(str(self.size()) + " alignment blocks were extracted from " + output_fname)

class SimpleBlastOutputConfig:
    header_str = "# BLAST"
    query_str = "# Query: "

class SimpleBlastOutput:
    match_map = dict()

    def LineIsQuery(self, line):
        return line[:len(SimpleBlastOutputConfig.query_str)] == SimpleBlastOutputConfig.query_str

    def ParseQueryName(self, line):
        return line[len(SimpleBlastOutputConfig.query_str):]

    def GetIdByName(self, name):
        return int(name.split("___")[1])

    def LineIsHit(self, line, query_name):
        return line[:len(query_name)] == query_name

    def GetHitSubject(self, line):
        splits = line.split('\t')
        return splits[1]

    def __init__(self, src_fname):
        self.match_map = dict()
        if not os.path.exists(src_fname):
            print("ERROR: BLAST output file " + src_fname + " was not found")
            sys.exit(1)
        src_fhandler = open(src_fname, 'r')

        query_name = ""
        for line in src_fhandler.readlines():
            line = line.strip()
            if self.LineIsQuery(line):
                query_name = self.ParseQueryName(line)
                self.match_map[self.GetIdByName(query_name)] = list()
            elif self.LineIsHit(line, query_name) and query_name != "" and len(self.match_map[self.GetIdByName(query_name)]) == 0:
                self.match_map[self.GetIdByName(query_name)].append(self.GetHitSubject(line))

        print(str(len(self.match_map)) + " blocks were extracted from " + src_fname)
