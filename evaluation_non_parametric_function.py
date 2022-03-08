import numpy as np
from scipy import stats
from scipy.integrate import quad
from patsy import bs
import statsmodels.api as sm

class evaluate_non_parametric():
    def __init__(self, x, y):
        self.x = x
        self.y = y
        # self.y_hat = y_hat
    
    def bs_estimate(self, knots, deg):
        bx = bs(x=self.x, knots=knots, degree=deg,
                include_intercept=True, lower_bound=0, upper_bound=1)
        model = sm.OLS(self.y, bx)
        result = model.fit()
        y_hat = bx @ result.params.reshape((-1,1))
        return y_hat
    
    def kernel_estimate(self, bandwidth, x0):
        y_hat = np.sum(self.y*stats.norm.pdf((x0-self.x)/bandwidth)) / np.sum(stats.norm.pdf((x0-self.x)/bandwidth))
        return y_hat
    
    def bs_ISE(self, x, knots, deg):
        bx = bs(x=x, knots=knots, degree=deg,
                include_intercept=True, lower_bound=0, upper_bound=1)
        
        m = lambda x: x*np.sin(20*x) + np.random.randn(1000)*0.2

        def m_hat(x):
            model = sm.OLS(self.y, bx)
            result = model.fit()
            m_hat_vector = bx @ result.params.reshape((-1,1))
            m_hat_vector = m_hat_vector.flatten()
            return m_hat_vector

        SE = lambda x: (m_hat(x) - m(x))**2

        ise, _ = quad(SE, 0, 1)
        return ise
