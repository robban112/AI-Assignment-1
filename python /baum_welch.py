import sys
import math
from matrix_utilities import Model, pretty_print, times_vectors, getcol, parse_model, matrix_to_string

## CREATE ALPHA MATRIX ##
def create_alpha_matrix2(model, obs_seq):
    a, b, pi, n = model
    last_t = len(obs_seq)
    alpha = [[0]*n for _ in range(0,last_t)]
    cs = [0]*last_t
    # compute alpha[0]
    c0 = 0
    for i in range(0,n):
        alpha[0][i] = pi[i]*b[i][obs_seq[0]]
        c0 += alpha[0][i]

    # scale the alpha[0]
    c0 = 1/c0
    cs[0] = c0
    for i in range(0,n):
        alpha[0][i] = c0*alpha[0][i]

    # compute alpha[t]
    for t in range(1,last_t):
        ct = 0
        for i in range(0,n):
            alpha[t][i] = 0
            for j in range(0,n):
                alpha[t][i] += alpha[t-1][j]*a[j][i]
            alpha[t][i] *= b[i][obs_seq[t]]
            ct += alpha[t][i]

        # scale alpha[t]
        ct = 1/ct
        cs[t] = ct
        for i in range(0,n):
            alpha[t][i] *= ct

    return (alpha, cs)

def create_alpha_matrix(model, obs_seq):
    a, b, pi, n = model
    alpha_0 = times_vectors(pi, getcol(b, obs_seq[0]))
    scaling_values = []
    #scale the alpha 0
    c0 = 1/sum(alpha_0)
    alpha_0 = [x*c0 for x in alpha_0]
    scaling_values.append(c0)
    ##

    alpha_matrix = [alpha_0]
    # loop through each of the observations
    for obs in obs_seq[1:]:
        prev_alpha_pass = alpha_matrix[-1]
        curr_alpha_pass = []
        # loop through states
        for j in range(0,n):
            obs_val = b[j][obs]
            if (obs_val == 0):
                 curr_alpha_pass.append(0)
                 continue
            trans_val = sum([prev_val * a[i][j] for (i,prev_val) in enumerate(prev_alpha_pass)])
            curr_alpha_pass.append(obs_val * trans_val)

        # scale alpha t (i)
        ct = 1/sum(curr_alpha_pass)
        scaling_values.append(ct)
        curr_alpha_pass = [x*ct for x in curr_alpha_pass]

        alpha_matrix.append(curr_alpha_pass)
    return (alpha_matrix, scaling_values)

## CREATE BETA MATRIX ##
def create_beta_matrix2(model, obs_seq, scaling_values):
    a, b, _, n = model
    last_t = len(obs_seq)
    beta = [[0]*n for _ in range(0,last_t)]
    # let beta[last_t-1] = 1, scaled by c_last_t-1
    for i in range(0,n):
        beta[last_t-1][i] = scaling_values[last_t-1]

    for t in range(last_t-2,-1,-1):
        for i in range(0,n):
            beta[t][i] = 0
            for j in range(0,n):
                beta[t][i] += a[i][j]*b[j][obs_seq[t+1]]*beta[t+1][j]
            # scale beta[t] with same scale factor as alpha[t]
            beta[t][i] *= scaling_values[t]
    return beta

def create_beta_matrix(model, obs_seq, scaling_values):
    a, b, pi, n = model
    beta_0 = [scaling_values[-1]]*n
    beta_matrix = [beta_0]
    last_t = len(obs_seq)
    # loop through each of the observations
    for t in range(last_t-2,-1,-1):
        prev_beta_pass = beta_matrix[0]
        curr_beta_pass = []
        # loop through states
        for i in range(0,n):
            val = 0
            for j, prev_val in enumerate(prev_beta_pass):
                val += b[j][obs_seq[t+1]] * a[i][j] * prev_val
            curr_beta_pass.append(val)

        #scale beta t (i)
        ct = scaling_values[t]
        curr_beta_pass = [x*ct for x in curr_beta_pass]
        beta_matrix.insert(0, curr_beta_pass)
    return beta_matrix

## CREATE DI_GAMMA AND GAMMA ##
def create_di_gamma_and_gamma_matrices(a, b, alpha, beta, obs_seq):
    n = len(a)
    last_t = len(obs_seq)
    di_gamma = [[[0]*n for _ in range(0,n)] for _ in range(0,last_t)]
    gamma =    [[0]*n for _ in range(0,last_t)]
    for t in range(0,last_t-1):
        denom = 0
        for i in range(0,n):
            for j in range(0,n):
                denom += alpha[t][i] * a[i][j] * b[j][obs_seq[t+1]] * beta[t+1][j]

        for i in range(0,n):
            gamma[t][i] = 0
            for j in range(0,n):
                di_gamma[t][i][j] = (alpha[t][i]*a[i][j]*b[j][obs_seq[t+1]]*beta[t+1][j])/denom
                gamma[t][i] += di_gamma[t][i][j]

    # Special case for gamma last t
    denom = 0
    for i in range(0,n):
        denom = denom + alpha[last_t-1][i]
    for i in range(0,n):
        gamma[last_t-1][i] = alpha[last_t-1][i]/denom
    return (gamma, di_gamma)

def create_di_gamma_matrix(a, b, alpha, beta, obs_seq):
    last_t = len(obs_seq)
    n = len(a)
    di_gamma_matrix = []
    for t in range(0,last_t-1):
        step = []
        denom = 0
        for i, j in [(i,j) for i in range(0,n) for j in range(0,n)]:
            denom += alpha[t][i]*a[i][j]*b[j][obs_seq[t+1]]*beta[t+1][j]

        for i in range(0,n):
            i_j_trans = []
            for j in range(0,n):
                val = alpha[t][i]*a[i][j]*b[j][obs_seq[t+1]]*beta[t+1][j]
                i_j_trans.append(val / denom)
            step.append(i_j_trans)

        di_gamma_matrix.append(step)
    return di_gamma_matrix

def create_gamma_matrix(di_gamma, alpha_matrix):
    # print("di_gamma: " + str(di_gamma))
    map_sum = lambda l: list(map(sum,l))
    last_sum = sum(alpha_matrix[-1])
    last_row = [x/last_sum for x in alpha_matrix[-1]]
    return [map_sum(step) for step in di_gamma] + [last_row]


## REESTIMATE PI ##
def reestimate_pi(gamma):
    return gamma[0]

#REESTIMATE A ##
def reestimate_a2(gamma, di_gamma):
    n = len(gamma[0])
    last_t = len(gamma)
    a = [[0]*n for _ in range(0,n)]
    for i in range(0,n):
        for j in range(0,n):
            numer = 0
            denom = 0
            for t in range(0,last_t-1):
                numer += di_gamma[t][i][j]
                denom += gamma[t][i]
            a[i][j] = numer/denom
    return a

def reestimate_a(gamma, di_gamma):
    print(gamma[:5])
    n = len(di_gamma[0])
    a = []
    for i in range(0, n):
        row = []
        for j in range(0, n):
            num = sum([step[i][j] for step in di_gamma])
            denom = sum([step[i] for step in gamma])
            if (num == 0): row.append(0)
            else: row.append(num / denom)
        a.append(row)
    return a

## REESTIMATE B ##
def reestimate_b2(gamma, obs_seq, m):
    n = len(gamma[0])
    last_t = len(obs_seq)
    b = [[0]*m for _ in range(0,n)]
    for i in range(0,n):
        for j in range(0,m):
            numer = 0
            denom = 0
            for t in range(0,last_t):
                if obs_seq[t] == j:
                    numer += gamma[t][i]
                denom += gamma[t][i]
            b[i][j] = numer/denom
    return b

def reestimate_b(gamma, obs_seq, m):
    n = len(gamma[0])
    b_matrix = []
    for j in range(0, n):
        b_row = []
        for k in range(0,m):
            num = 0
            denom = 0
            for t, obs in enumerate(obs_seq):
                denom += gamma[t][j]
                if (obs == k):
                    num += gamma[t][j]
            #if denom == 0: b_row.append(0)
            b_row.append(num / denom)
        b_matrix.append(b_row)
    return b_matrix

## REESTIMATE MODEL ##
def reestimate_model(model, obs_seq, alpha, beta):
    gamma, di_gamma = create_di_gamma_and_gamma_matrices(model.a, model.b, alpha, beta, obs_seq)
    # di_gamma = create_di_gamma_matrix2(model.a, model.b, alpha, beta, obs_seq)
    # gamma = create_gamma_matrix(di_gamma, alpha)
    new_pi = reestimate_pi(gamma)
    new_a = reestimate_a2(gamma, di_gamma)
    new_b = reestimate_b2(gamma, obs_seq, len(model.b[0]))
    return Model(a=new_a, b=new_b, pi=new_pi, n=len(new_a))

## BAUM WELCH ##
def baum_welch2(model, obs_seq):
    maxIters = 30
    iters = 0
    oldLogProb = -10000
    logProb = -9999
    while iters < maxIters and logProb > oldLogProb:
        oldLogProb = logProb
        alpha, cs = create_alpha_matrix2(model, obs_seq)
        beta = create_beta_matrix2(model, obs_seq, cs)
        model = reestimate_model(model, obs_seq, alpha, beta)
        logProb = compute_log_prob(cs)
        iters += 1
    return model

def baum_welch(model, obs_seq):
    prev_alpha, prev_scaling_values = create_alpha_matrix2(model, obs_seq)
    prev_beta = create_beta_matrix2(model, obs_seq, prev_scaling_values)
    n = 0
    maxIter = 80
    logProb, oldLogProb = 1, 0
    while n < maxIter and logProb > oldLogProb:
        reestimated_model = reestimate_model(model, obs_seq, prev_alpha, prev_beta)
        new_alpha, new_scaling_values = create_alpha_matrix2(reestimated_model, obs_seq)
        new_beta = create_beta_matrix2(reestimated_model, obs_seq, new_scaling_values)
        oldLogProb = compute_log_prob(prev_scaling_values)
        logProb = compute_log_prob(new_scaling_values)
        #print("logProb: " + str(logProb))
        #print("newLog: " + str(oldLogProb))
        prev_alpha = new_alpha
        prev_beta = new_beta
        prev_scaling_values = new_scaling_values
        model = reestimated_model
        n += 1
        print("n: " + str(n))
    return model

## LOG PROB ##
def compute_log_prob(scaling_values):
    logProb = sum(map(math.log, scaling_values))
    return -logProb

## TEST ##
def test():
    a = [[0.3, 0.7],[0.9, 0.1]]
    b = [[0.5,0.5],[0.5,0.5]]
    pi = [1,0]
    n = 2
    obs_seq = [0, 0, 1, 1]
    model = Model(a, b, pi, n)
    print(" --- A --- ")
    pretty_print(a)
    print(" --- B --- ")
    pretty_print(b)
    print(" --- pi --- ")
    print(pi)
    print(baum_welch(model, obs_seq))

def parse_observation_seq():
    line = sys.stdin.readline()
    seq = line.strip().split()[1:]
    return list(map(int, seq))

def print_output(model):
    print(matrix_to_string(model.a))
    print(matrix_to_string(model.b))

## TOP LEVEL CODE ##

model = parse_model()
obs_seq = parse_observation_seq()

baum_model = baum_welch2(model, obs_seq)
print_output(baum_model)
