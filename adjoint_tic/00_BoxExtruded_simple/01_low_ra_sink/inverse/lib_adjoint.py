class myReducedFunctional(ReducedFunctional):
    def __init__(self, functional, controls, **kwargs):
        super().__init__(functional, control, **kwargs)
        self.fwd_cntr = 0 
        self.adj_cntr = 0 

    def __call__(self, values):
        pre_time = time.time()
        val = super().__call__(values)
        self.fwd_cntr += 1
        log(f"\tFWD # {self.fwd_cntr} call took {time.time() - pre_time}")
        return val 

    def derivative(self):
        pre_time = time.time()
        deriv = super().derivative(options={})
        self.adj_cntr += 1
        log(f"\tADJ # {self.adj_cntr} call took {time.time() - pre_time}")
        return deriv

