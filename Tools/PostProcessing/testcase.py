class TestCase:

    def __init__( self, name ):
        self.name = name
        if ( name == 'pml_x_psatd' ):
           self.energy_ref = 9.1301289517e-08
           self.reflec_ref = 1.3806831258153887e-06
        if ( name == 'pml_x_yee' ):
           self.energy_ref = 9.1301289517e-08
           self.reflec_ref = 5.683000058954201e-07
        if ( name == 'pml_x_ckc' ):
           self.energy_ref = 9.1301289517e-08
           self.reflec_ref = 1.8015e-06
