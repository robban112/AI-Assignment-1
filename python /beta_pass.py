import sys
import math
from matrix_utilities import Model, pretty_print, times_vectors, getcol, parse_model

def create_alpha_matrix(model, obs_seq):
    a, b, pi, n = model
    alpha_0 = times_vectors(pi, getcol(b, obs_seq[0]))
    scaling_values = []
    #scale the alpha 0
    c0 = 1/math.fsum(alpha_0)
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
            trans_val = math.fsum([prev_val * a[i][j] for (i,prev_val) in enumerate(prev_alpha_pass)])
            curr_alpha_pass.append(obs_val * trans_val)

        # scale alpha t (i)
        ct = 1/math.fsum(curr_alpha_pass)
        scaling_values.append(ct)
        curr_alpha_pass = [x*ct for x in curr_alpha_pass]

        alpha_matrix.append(curr_alpha_pass)
    return (alpha_matrix, scaling_values)

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

def di_gamma(t, i, j, a, b, alpha, beta, obs_seq):
    num = alpha[t][i]*a[i][j]*b[j][obs_seq[t+1]]*beta[t+1][j]
    denom = math.fsum(alpha[len(obs_seq)-1])
    #if (denom == 0):
    #   return 0
    return num / denom

def create_di_gamma_matrix2(a, b, alpha, beta, obs_seq):
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

def create_di_gamma_matrix(a, b, alpha, beta, obs_seq):
    di_gamma_matrix = []
    t = len(obs_seq)
    for t in range(0,t-1):
        t_trans = []
        for i in range(0,len(a)):
            i_j_trans = []
            for j in range(0,len(a)):
                i_j_trans.append(di_gamma(t, i, j, a, b, alpha, beta, obs_seq))
            t_trans.append(i_j_trans)
        di_gamma_matrix.append(t_trans)
    return di_gamma_matrix

# def gamma(di_gamma_matrix, t, i):
#     return math.fsum(di_gamma_matrix[t][i])

def create_gamma_matrix(di_gamma, alpha_matrix):
    # print("di_gamma: " + str(di_gamma))
    map_sum = lambda l: list(map(math.fsum,l))
    last_sum = math.fsum(alpha_matrix[-1])
    last_row = [x/last_sum for x in alpha_matrix[-1]]
    return [map_sum(step) for step in di_gamma] + [last_row]

def reestimate_pi(gamma):
    return gamma[0]

def reestimate_a(gamma, di_gamma):
    n = len(di_gamma[0])
    a = []
    for i in range(0, n):
        row = []
        for j in range(0, n):
            num = math.fsum([step[i][j] for step in di_gamma])
            denom = math.fsum([step[i] for step in gamma[:-1]])
            #print("di: " + str(di_gamma))
            #print("ga: " + str(gamma))
            #print("len di: " + str(len(di_gamma)))
            #print("len gam: " + str(len(gamma)))
            #if (denom == 0): row.append(0)
            #print("i,j: " + str(i) + ", " + str(j))
            #print("num: " + str(num) + ", denom: " + str(denom))
            if (num == 0): row.append(0)
            else: row.append(num / denom)
        a.append(row)
    return a

def reestimate_b(gamma, obs_seq, m):
    n = len(gamma[0])
    b_matrix = []
    for j in range(0, n):
        b_row = []
        for k in range(0,m):
            num = 0
            denom = 0
            for idx, obs in enumerate(obs_seq):
                denom += gamma[idx][j]
                if (obs == k):
                    num += gamma[idx][j]
            #if denom == 0: b_row.append(0)
            b_row.append(num / denom)
        b_matrix.append(b_row)
    return b_matrix

def reestimate_model(model, obs_seq, alpha, beta):
    di_gamma = create_di_gamma_matrix2(model.a, model.b, alpha, beta, obs_seq)
    gamma = create_gamma_matrix(di_gamma, alpha)
    new_pi = reestimate_pi(gamma)
    new_a = reestimate_a(gamma, di_gamma)
    new_b = reestimate_b(gamma, obs_seq, len(model.b[0]))
    # print(" --- alpha ---")
    # print(alpha)
    # print(" --- beta ---")
    # pretty_print(beta)
    # print(" --- di gamma ---")
    # print(di_gamma)
    # print(" --- reestimated pi --- ")
    # print(new_pi)
    # print(" --- reestimated a --- ")
    # pretty_print(new_a)
    # print(" --- reestimated b --- ")
    # pretty_print(new_b)
    return Model(a=new_a, b=new_b, pi=new_pi, n=len(new_a))

def baum_welch(model, obs_seq):
    prev_alpha, prev_scaling_values = create_alpha_matrix(model, obs_seq)
    prev_beta = create_beta_matrix(model, obs_seq, prev_scaling_values)
    n = 0
    maxIter = 200
    logProb, oldLogProb = 1, 0
    while maxIter > n and logProb > oldLogProb:
        reestimated_model = reestimate_model(model, obs_seq, prev_alpha, prev_beta)
        new_alpha, new_scaling_values = create_alpha_matrix(reestimated_model, obs_seq)
        new_beta = create_beta_matrix(reestimated_model, obs_seq, new_scaling_values)
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

def compute_log_prob(scaling_values):
    logProb = math.fsum(map(math.log, scaling_values))
    return -logProb

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

model = parse_model()
obs_seq = parse_observation_seq()

print(baum_welch(model, obs_seq))
