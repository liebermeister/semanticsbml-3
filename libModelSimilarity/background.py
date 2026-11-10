import numpy
import os, os.path, random, scipy, pickle, libsbml, math, scipy.misc, datetime
from . import compareModels
from functools import cmp_to_key

random_models_filename = os.path.join( os.path.expanduser('~'), '.semanticSBML', 'random_biomodels' )
biomodels_probabilities_filename = os.path.join( os.path.expanduser('~'), '.semanticSBML', 'biomodels_probabilities' )

class Background(compareModels.CompareModels):
    def __init__(self, model, directory):
        compareModels.CompareModels.__init__(self)
        if type(model)==type('str'):
            self.model = compareModels.ComparableModel(self, filename=model)
        elif type(model) == libsbml.SBMLDocument:
            self.model = compareModels.ComparableModel(self, sbmldoc=model)
        else:
            self.model = compareModels.ComparableModel(self, model_elements_annotations=model)
        self.model_vector = numpy.matrix(self.model._model_elements_annotation_to_element_vector(self.model.mea)).transpose()
        self.directory = directory
        self.freqs = None
        self.scores = None
        self.beta_density = None
        self.simple_binomial = None
    def _read_dir(self, directory):
        filenames = os.listdir(directory)
        filenames.sort()
        complete_filenames = []
        for filename in filenames:
            if filename.startswith('BIOMD') and filename.endswith('.xml'):
                complete_filenames.append( os.path.join(directory, filename) )
        return complete_filenames
    def _filename_to_vector(self, filename):
        cm = compareModels.ComparableModel(self, filename=filename)
        vector = cm._model_elements_annotation_to_element_vector(cm.mea)
        return vector
    def _vector_to_problist(self, vector):
        problist = []
        for index in range(len(vector)):
            if vector[index] > .000000001:
                problist.append(self.freqs[index])
        problist.sort()
        return problist
    def _get_frequencies(self, directory):
        filenames = self._read_dir(directory)
        vectors = [self._filename_to_vector(f) for f in filenames]
        freqs = [1.]*len(self.annotation_keys)
        for i in range(len(vectors)):
            for j in range(len(vectors[i])):
                freqs[j] = freqs[j] + vectors[i][j]
        for j in range(len(freqs)):
            freqs[j] = freqs[j] / float(len(vectors)+1)
        return freqs
    def get_random_vector(self):
        return numpy.matrix([1. if random.random()<prob else 0. for prob in self.freqs]).transpose()
    def write_random_matrix(self,size):
        matrix = [self.get_random_vector() for i in range(size)]
        #f=open(random_models_filename,'w')
        #pickle.dump(matrix,f)
        #f.close()
        return matrix
    def get_random_matrix(self,size):
        if False and os.path.isfile(random_models_filename):
            f=open(random_models_filename,'r')
            matrix = pickle.load(f)
            f.close()
        else:
            matrix = self.write_random_matrix(size)
        return matrix
    def _load_freqs(self):
        if os.path.isfile(biomodels_probabilities_filename):
            f = open(biomodels_probabilities_filename,'r')
            self.freqs = pickle.load(f)
            f.close()
        else:
            self.freqs = self._get_frequencies(self.directory)
            f = open(biomodels_probabilities_filename,'w')
            pickle.dump(self.freqs,f)
            f.close()
    def get_scores(self,size):
        if self.freqs == None:
            self._load_freqs()
        scores = []
        for i in range(size):
            scores.append(self.similarity(self.model_vector, self.get_random_vector()))
        self.scores = scores
        return scores
    def get_p_scores(self,sim):
        if self.scores == None:
            self.get_scores(10000)
        return float(len([x for x in self.scores if x > sim])+1)/float(len(self.scores)+2)
    def get_p_beta(self,sim):
        if self.scores == None:
            self.get_scores(10000)
        if self.beta_density == None:
            self.get_beta_dist(self.scores)
        [betax, betay] = self.beta_density
        list_of_smaller = [x for x in betax if x < sim]
        if len(list_of_smaller) > 0:
            close_to_sim = max(list_of_smaller)
        else:
            close_to_sim = min(betax)
        return betay[betax.index(close_to_sim)]
    def get_p_binomial(self, sim, n):
        # does not work yet  m(
        if self.freqs == None:
            self._load_freqs()
        a = len(self.freqs)
        p = sum(self.freqs)/float(a)
        m = self.model_vector.sum()
        xmax = int(math.floor(sim*math.sqrt(m)*math.sqrt(n)))
        if self.simple_binomial is None:
            self.simple_binomial = [(1.-p)**a]
        k=1
        while len(self.simple_binomial) <= xmax:
            self.simple_binomial.append( self.simple_binomial[-1] + scipy.special.comb(a, k) * p**k * (1.-p)**(a-k) )
            k = k + 1
        return 1.-self.simple_binomial[xmax]
    def precalc_p_dynamic(self,max_overlap):
        import scipy.signal
        bound = max(int(math.floor(max_overlap)),2)
        m = int(self.model_vector.sum())
        this_freqs = [self.freqs[i] for i in range(len(self.freqs)) if self.model_vector[i][0] > 0.]
        # build the dynamic programming matrix
        self.A = numpy.array([1.])
        for i in range(m):
            self.A = scipy.signal.convolve(self.A, numpy.array([1.-this_freqs[i], this_freqs[i]]))[:bound]
    def get_p_dynamic(self, overlap):
        import scipy.signal
        bound = int(math.floor(overlap))
        return 1.-sum(self.A[:bound])
    def get_p_wolf_exp(self, sim):
        m = int(self.model_vector.sum())
        start = datetime.datetime.now()
        matrix = numpy.eye(m+1)
        this_freqs = [self.freqs[i] for i in range(len(self.freqs)) if self.model_vector[i][0] > 0.]
        # build the dynamic programming matrix
        for i in range(m):
            # fill a row in the matrix
            freq = this_freqs[i]
            invf = 1.-freq
            matrix[i+1][i+1] = 0.
            for j in range(i+1):
                matrix[i+1][j] = matrix[i+1][j] + invf*matrix[i][j]
                matrix[i+1][j+1] = matrix[i+1][j+1] + freq*matrix[i][j]
        A = matrix[m]
        print("A",datetime.datetime.now()-start)
        start = datetime.datetime.now()
        other_freqs = [self.freqs[i] for i in range(len(self.freqs)) if self.model_vector[i][0] == 0.]
        prob_and_count = [[x, other_freqs.count(x)] for x in set(other_freqs)]
        prob_and_count.sort(key=cmp_to_key(lambda x,y:x[1]-y[1]))
        B = [1.]
        for [p, a] in prob_and_count:
            binomial = [scipy.special.comb(a, k) * p**k * (1.-p)**(a-k) for k in range(a+1)]
            new_B = [0.]*(len(B)+len(binomial)-1)
            for i in range(len(B)):
                for j in range(len(binomial)):
                    new_B[i+j] = new_B[i+j] + B[i]*binomial[j]
            B = new_B
        print("B",datetime.datetime.now()-start)
        start = datetime.datetime.now()
        prob_sum = 0.
        for a in range(1,len(A)):
            for b in range(len(B)):
                score = float(a)/math.sqrt(a+b)
                if score < sim:
                    break
                prob_sum += A[a]*B[b]
        print("sum",datetime.datetime.now()-start)
        return prob_sum
    def get_p_wolf(self, sim):
        m = int(self.model_vector.sum())
        start = datetime.datetime.now()
        import scipy
        matrix = numpy.eye(m+1)
        this_freqs = [self.freqs[i] for i in range(len(self.freqs)) if self.model_vector[i][0] > 0.]
        # build the dynamic programming matrix
        for i in range(m):
            # fill a row in the matrix
            freq = this_freqs[i]
            invf = 1.-freq
            matrix[i+1][i+1] = 0.
            for j in range(i+1):
                matrix[i+1][j] = matrix[i+1][j] + invf*matrix[i][j]
                matrix[i+1][j+1] = matrix[i+1][j+1] + freq*matrix[i][j]
        A = matrix[m]
        print("A",datetime.datetime.now()-start)
        start = datetime.datetime.now()
        # build B from normal distribution
        mu = 0.
        sigma2 = 0.
        other_freqs = [self.freqs[i] for i in range(len(self.freqs)) if self.model_vector[i][0] == 0.]
        prob_and_count = [[x, other_freqs.count(x)] for x in set(other_freqs)]
        for [p, n] in prob_and_count:
            mu += n*p
            sigma2 += n*p*(1.-p)
        import scipy.stats.distributions
        normal = scipy.stats.distributions.norm(loc=mu, scale=math.sqrt(sigma2))
        cdfs = [normal.cdf(i-.5) for i in range(len(other_freqs)+2)]
        B = [cdfs[i+1] - cdfs[i] for i in range(len(other_freqs)+1)]
        print("B",datetime.datetime.now()-start)
        start = datetime.datetime.now()
        prob_sum = 0.
        for a in range(1,len(A)):
            for b in range(len(B)):
                score = float(a)/math.sqrt(a+b)
                if score < sim:
                    break
                prob_sum += A[a]*B[b]
        print("sum",datetime.datetime.now()-start)
        return prob_sum
    def get_alpha_beta(self, scores):
        mean = sum(scores)/float(len(scores))
        variance = sum([(x-mean)**2 for x in scores])/float(len(scores))
        alpha = mean*( (mean*(1.-mean))/variance - 1. )
        beta = (1.-mean)*( (mean*(1.-mean))/variance - 1. )
        return [alpha, beta, mean, variance]
    def get_beta_dist(self, scores, bins=30):
        mean = sum(scores)/float(len(scores))
        variance = sum([(x-mean)**2 for x in scores])/float(len(scores))
        alpha = mean*( (mean*(1.-mean))/variance - 1. )
        beta = (1.-mean)*( (mean*(1.-mean))/variance - 1. )
        def beta_func(norm,a,b,x):
            return 1./norm * x**(a-1.) * (1.-x)**(b-1.)
        num_ints = 1000
        len_int = 1./num_ints
        int = 0.
        for i in range(num_ints):
            x = i*len_int + .5*len_int
            int += len_int * beta_func(1., alpha, beta, x)
        outx=[]
        outy=[]
        cum_outx=[]
        cum_outy=[]
        cumulative = 0.
        inth = max(scores)/float(bins)*len(scores)
        for i in range(num_ints):
            x = i*len_int + .5*len_int
            y = beta_func(int/inth, alpha, beta, x)
            outx.append(x)
            outy.append(y)
        norm = sum(outy)
        for i in range(num_ints):
            cumulative = cumulative + outy[num_ints-i-1]/norm
            cum_outx.append(outx[num_ints-i-1])
            cum_outy.append(cumulative)
        self.beta_density = [cum_outx, cum_outy]
        return [outx,outy]
def main():
    number = 1000
    bins = 30
    import datetime
    from pylab import hist
    import matplotlib

    path = '/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/'

    for model_number in range(1,232):
        start = datetime.datetime.now()
    
        fig=matplotlib.pyplot.figure(frameon=False)
        ax = fig.add_subplot(111)

        filename = path+'BIOMD'+str(model_number).rjust(10,'0')+'.xml'
        bg = Background(filename,path)
        scores = bg.get_scores(number)

        try:
            [alpha, beta, mean, variance] = bg.get_alpha_beta(scores)
            problist = bg._vector_to_problist(bg._filename_to_vector(filename))
            len_prob = len(problist)
            sum_prob = sum(problist)
            mean_prob = sum_prob/len_prob

            [betax,betay] = bg.get_beta_dist(scores,bins)
            f=file('blabla.csv','a')
            f.write(','.join([str(x) for x in [alpha,beta,len_prob,sum_prob,mean_prob,1.-(sum(betay[:100])/sum(betay)),1.-(sum(betay[:250])/sum(betay)),1.-(sum(betay[:500])/sum(betay)),mean,variance]])+'\n')
            f.close()
            print(filename)
        except ZeroDivisionError:
            continue
        continue
        (y,x,p)=hist(scores,bins)
        ax.plot(betax,betay)
        fig.savefig('/home/schulzma/Desktop/p/sim'+str(model_number)+'.png')
        print('took',datetime.datetime.now()-start)

if __name__ == "__main__":
    #main()
    path = '/home/schulzma/semanticsbml/trunk/libModelSimilarity/curated/'
    filename = path+'BIOMD0000000061.xml'
    bg = Background(filename,path)
    #for sim in [.01,.02,.05,.1,.12,.15,.2,.25,.5]:
    #    print sim, bg.get_p_scores(sim), bg.get_p_beta(sim)
    bg.get_scores(998)
    import datetime
    start = datetime.datetime.now()
    bg.precal_p_dynamic(60)
    for sim in reversed(list(range(60))):
        print(sim, bg.get_p_dynamic(sim))
    print(datetime.datetime.now()-start)
