import time
from biopylib.sequence import SQ

# Heuristic Method for deterministic motif identification

class hdMotif:
    
    # Constructor
    def __init__(self, size = 3, seqs = None):
        self.motif_size = size
        if (seqs is not None):
            self.seqs = seqs
        else:
            self.seqs = []
            
    # Size of sequence
    def seq_size(self,seq_id):
        return len(self.seqs[seq_id])
            
    # Consensus (heuristic)   
    def heuristic_consensus(self):
        res = [0]*len(self.seqs) 
        max_score = -1
        partial = [0,0]
        for i in range(self.seq_size(0)-self.motif_size):
            for j in range(self.seq_size(1)-self.motif_size):
                partial[0] = i
                partial[1] = j
                sc = self.score(partial)
                if(sc > max_score):
                    max_score = sc
                    res[0] = i; res[1] = j
        for k in range(2, len(self.seqs)):
            partial = [0]*(k+1)
            for j in range(k):
                partial[j] = res[j]
            max_score = -1
            for i in range(self.seq_size(k)-self.motif_size):
                partial[k] = i
                sc = self.score(partial)
                if(sc > max_score):
                    max_score = sc
                    res[k] = i
        return res    
    
# Exhaustive Search Motif Discovery
    
class edMotif:
    
    # Constructor
    def __init__(self, size = 3, seqs = None):
        self.motif_size = size
        if (seqs is not None):
            self.seqs = seqs
        else:
            self.seqs = []
            
    # Size of sequence
    def seq_size(self,seq_id):
        return len(self.seqs[seq_id])
    
    def next_solution (self, s):
        next_sol= [0]*len(s)
        pos = len(s) - 1     
        while(pos >=0 and s[pos] == self.seq_size(pos) - self.motif_size):
            pos -= 1
        if (pos < 0): 
            next_sol = None
        else:
            for i in range(pos): 
                next_sol[i] = s[i]
            next_sol[pos] = s[pos]+1
            for i in range(pos+1, len(s)):
                next_sol[i] = 0
        return next_sol

    def exhaustive_search(self):
        best_score = -1
        res = []
        s = [0]* len(self.seqs)
        while (s!= None):
            sc = self.score(s)
            if (sc > best_score):
                best_score = sc
                res = s
            s = self.next_solution(s)
            
        return res
            
# Branch and Bound Motif Discovery Class 

class bdMotif:

    # Constructor
    def __init__(self, size = 3, seqs = None):
        self.motif_size = size
        if (seqs is not None):
            self.seqs = seqs
        else:
            self.seqs = [] 
            
    # Size of sequence
    def seq_size(self,seq_id):
        return len(self.seqs[seq_id])
            
    def next_vertex(self, s):
        res =  []
        if len(s) < len(self.seqs): # internal node -> down one level
            for i in range(len(s)): 
                res.append(s[i])
            res.append(0)
        else: 
            pos = len(s)-1 
            while pos >=0 and s[pos] == self.seq_size(pos) - self.motif_size:
                pos -= 1
            if pos < 0: res = None # last solution
            else:
                for i in range(pos): res.append(s[i])
                res.append(s[pos]+1)
        return res
    
    def bypass(self, s):
        res =  []
        pos = len(s)-1
        while(pos >=0 and s[pos] == self.seq_size(pos) - self.motif_size):
            pos -= 1
        if(pos < 0):
            res = None 
        else:
            for i in range(pos): res.append(s[i])
            res.append(s[pos]+1)
        return res
     
    def branch_and_bound(self):
        best_score = -1
        best_motif = None
        size = len(self.seqs)
        s = [0]*size
        while s != None:
            if len(s) < size:
                # estimate the bound for current internal node
                # test if the best score can be reached
                optimum_score = self.score(s) + (size-len(s)) * self.motif_size
                if optimum_score < best_score: s = self.bypass(s)
                else: s = self.next_vertex(s)
            else:
                # test if current leaf is a better solution
                sc = self.score(s)
                if sc > best_score:
                    best_score = sc
                    best_motif = s
                s = self.next_vertex(s)
                
        return best_motif
        

''' 

Deterministic Motif Discovery Main Class  
- Heuristic Approach
- Exhaustive Approach (brute force)
- Branch and Bound Approach

'''

class dMotif(hdMotif,edMotif,bdMotif):

    # Constructor
    def __init__(self, size = 8, method='heurisic', seqs = None):
        self.motif_size = size
        if (seqs is not None):
            self.seqs = seqs
            self.alphabet = seqs[0].count.abc()
        else:
            self.seqs = []
            self.alphabet = "ACGT" # default: DNA
            
        self.method = method  # approach to dmd
        print(f'[note] using {self.method} approach')
        
        # Get access to class methods
        if(self.method == 'heuristic'):
            
            hdMotif.__init__(self,
                             size=self.motif_size,
                             seqs=self.seqs)
            
        elif(self.method == 'exhaustive'):
            
            edMotif.__init__(self,
                             size=self.motif_size,
                             seqs=self.seqs)
            
        elif(self.method == 'branch'):
            
            bdMotif.__init__(self,
                             size=self.motif_size,
                             seqs=self.seqs)
            
    # Special Methods
        
    def __len__ (self):
        return len(self.seqs)
    
    def __getitem__(self, n):
        return self.seqs[n]
    
    def seq_size(self,i):
        return len(self.seqs[i])
    
    ''' 
    
    Simple File Reader 
    sequences separated by new lines 
    
    '''
    
    def read_aln(self, fic, t):
        for s in open(fic, "r"):
            self.seqs.append(SQ(s.strip().upper(),t))
        self.alphabet = self.seqs[0].count.abc() # need to set
        
    def create_motif_from_indexes(self, indexes):
        res = [[0]*self.motif_size for i in range(len(self.alphabet))]
        for i,ind in enumerate(indexes):
            subseq = self.seqs[i][ind:(ind+self.motif_size)]
            for i in range(self.motif_size):
                for k in range(len(self.alphabet)):
                    if subseq[i] == self.alphabet[k]:
                        res[k][i] = res[k][i] + 1
        return res    

    def score(self,s):
        score = 0
        mat = self.create_motif_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score += maxcol
        return score
   
    def score_multiplicative(self, s):
        score = 1.0
        mat = self.create_motif_from_indexes(s)
        for j in range(len(mat[0])):
            maxcol = mat[0][j]
            for  i in range(1, len(mat)):
                if mat[i][j] > maxcol: 
                    maxcol = mat[i][j]
            score *= maxcol
        return score
    
    # Get the solution 
    
    def get(self):
        
        t0 = time.time()
        if(self.method == 'heuristic'):
            res = self.heuristic_consensus()
        elif(self.method == 'exhaustive'):
            res = self.exhaustive_search()
        elif(self.method == 'branch'):
            res = self.branch_and_bound()
            
        t1 = time.time()
        print(f'[note] finished {t1-t0:.3f} s')
            
        return (res,self.score(res))
