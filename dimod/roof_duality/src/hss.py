from dimod import DiscreteQuadraticModel
from itertools import combinations
import numpy as np


def get_problem_dimod(nv, nc, p, seed):
    np.random.seed(seed)
    dqm = DiscreteQuadraticModel()
    for u in range(nv):
        u = dqm.add_variable(nc)
        dqm.set_linear(u, np.random.normal(0.0, 3.0, size=nc))
    edges = [(i, j) for i, j in combinations(range(nv), r=2) if np.random.binomial(1, p)]
    print(len(edges))
    for u, v in edges:
        dqm.set_quadratic(u, v, np.random.normal(0.0, 3.0, size=[nc, nc]))
    return dqm


if __name__ == '__main__':
    from dqm_solver import HSSDQMSampler
    sampler = HSSDQMSampler()
    dqm = get_problem_dimod(100, 4, p=1.0, seed=1)
    sampler.sample(dqm, time_limit=160)
