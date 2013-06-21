package quantscale.analytic

class NIGModelSpec(val alpha: Double, val beta: Double, val delta: Double, val mu: Double) {
}

object NIGVanilla {
    def makeNIGModelSpecFromSABR(sabr : SABRModelSpec, forward: Double, tte: Double) : NIGModelSpec = {
      val s0 = sabr.alpha
      val rho = sabr.rho
      val nu = sabr.nu
      
      val A = (1+4*rho*rho)/5
      val B = -2*rho*rho
      val x = math.exp(nu*nu*tte)
      val kurtosis = A*(x*x*x*x+2*x*x*x+3*x*x+4*x+5)+2*B*(x+2)
      val variance = s0*s0*(x-1)/(nu*nu) 
      val deltaGamma = 3*(1+4*rho*rho)/(kurtosis-3)
      val alpha = math.sqrt(deltaGamma/variance)/(1-rho*rho)
      val beta = rho*alpha
      val gammaSq = (alpha*alpha-beta*beta)
      val mu = forward-deltaGamma*beta/gammaSq
      val delta = mu/math.sqrt(gammaSq)
      return new NIGModelSpec(alpha, beta, delta, mu)
    }
  
}