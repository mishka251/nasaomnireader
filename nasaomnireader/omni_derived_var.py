class omni_derived_var(object):
    """
    A variable that can be defined as a function
    of CDF variables, that we would like to be accessible via the same __getitem__
    interface as the variables contained in the CDF
    """

    def __init__(self, oi):
        self.oi = oi  # associated omni_interval
        self.varvals = None  # To hold computed so we don't recompute without needing to
        self.attrs = dict()