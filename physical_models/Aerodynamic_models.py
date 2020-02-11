import math
import random
import config

def aerodynamicscoefficient_constant(S, a):

    CL_body = 0
    CD_body = 2
    return CL_body, CD_body


def aerodynamicscoefficient_fM (T, S, a, body, vel_pp_mag, rho, gamma,aoa,montecarlo):
    alpha = aoa
    a_f = a.accomodation_factor
    Tw = T
    # Aerodynamics
    #if gamma > 0.2:
    #   alpha = math.pi/2

    def pressure(S,aoa,rho_inf,velocity,sigma):
        p = (rho_inf*velocity**2)/(2*S**2)*((((2-sigma)/math.pi**0.5)*S*math.sin(aoa)+((T/Tw)**0.5)*sigma/2)*math.exp(-(S*math.sin(aoa))**2)+((2-sigma)*((S*math.sin(aoa))**2+1/2)+sigma/2*math.pi**0.5*(S*math.sin(aoa)))*(1+math.erf(S*math.sin(aoa))))
        return p

    def tau(S,aoa,rho_inf,velocity,sigma):
        t = ((sigma * math.cos(aoa) * rho_inf*velocity**2) / (math.pi ** 0.5 * 2* S)) * (
                math.exp(-(S * math.sin(aoa)) ** 2) + math.pi ** 0.5 * (S * math.sin(aoa)) * (
                    1 + math.erf(S * math.sin(aoa))))
        return t

    # N = pressure(S, alpha, rho, vel_pp_mag, a_f )*body.Area_SA + pressure(S, math.pi/2, rho, vel_pp_mag, a_f)*body.Area_SC
    # A = tau(S, alpha, rho, vel_pp_mag, a_f )*body.Area_SA + tau(S, math.pi/2, rho, vel_pp_mag, a_f)*body.Area_SC

    # q = 1 / 2 * rho * vel_pp_mag ** 2
    # CN = N/(q*(body.Area_SA + body.Area_SC))
    # CA = A/(q*(body.Area_SA + body.Area_SC))
    def normalcoefficient(S,aoa,sigma):
        CN = 1 / (S ** 2) * ((((2 - sigma) / math.pi ** 0.5) * S * math.sin(aoa) + sigma / 2) * math.exp(
            -(S * math.sin(aoa)) ** 2) + (
                                         (2 - sigma) * ((S * math.sin(aoa)) ** 2 + 1 / 2) + sigma / 2 * math.pi ** 0.5 * (
                                             S * math.sin(aoa))) * (1 + math.erf(S * math.sin(aoa))))
        return CN

    def axialcoefficient(S,aoa,sigma):
        CA = ((sigma * math.cos(aoa)) / (math.pi ** 0.5 * S)) * (
                    math.exp(-(S * math.sin(aoa)) ** 2) + math.pi ** 0.5 * (S * math.sin(aoa)) * (
                        1 + math.erf(S * math.sin(aoa))))
        return CA


    #CN = 1/(S**2)*((((2-a_f)/math.pi**0.5)*S*math.sin(alpha)+a_f/2)*math.exp(-(S*math.sin(alpha))**2)+((2-a_f)*((S*math.sin(alpha))**2+1/2)+a_f/2*math.pi**0.5*(S*math.sin(alpha)))*(1+math.erf(S*math.sin(alpha))))
    #CA = ((a_f*math.cos(alpha))/(math.pi**0.5*S))*(math.exp(-(S*math.sin(alpha))**2)+math.pi**0.5*(S*math.sin(alpha))*(1+math.erf(S*math.sin(alpha))))
    ## Solar Panels:
    CN_sa = normalcoefficient(S,alpha,a_f)
    CA_sa = axialcoefficient(S,alpha,a_f)
    CL_sa = CN_sa*math.cos(alpha)-CA_sa*math.sin(alpha)
    CD_sa = CA_sa*math.cos(alpha)+CN_sa*math.sin(alpha)

    ## Spacecraft
    CN_sc = normalcoefficient(S,math.pi/2,a_f)
    CA_sc = axialcoefficient(S,math.pi/2,a_f)
    CL_sc = CN_sc*math.cos(math.pi/2)-CA_sc*math.sin(math.pi/2)
    CD_sc = CA_sc*math.cos(math.pi/2)+CN_sc*math.sin(math.pi/2)

    CD_body = (CD_sa*body.Area_SA + CD_sc*body.Area_SC)/(body.Area_SA+body.Area_SC)
    CL_body = (CL_sa*body.Area_SA + CL_sc*body.Area_SC)/(body.Area_SA+body.Area_SC)

    if montecarlo == True:
        random.seed(int(config.index_MonteCarlo))
        r = random.uniform(-CD_body*.1, CD_body*.1) # uncertanties 10%
        s = random.uniform(-CL_body*.1, CL_body*.1) # uncertanties 10%
        CD_body = CD_body+r
        CL_body = CL_body+s


    return CL_body, CD_body

## this is not done correctly
def aerodynamicscoefficient_noballisticflight(NoseRadius, BaseRadius , delta , alpha):
    theta = math.pi / 2 - delta

    # Probe
    CD_nose = 1 - (math.cos(theta))**4
    CL_nose = 0

    CN_conicfrostrum = (1-(NoseRadius/BaseRadius)**2 * (math.cos(delta)**2)) * (math.cos(delta))**2 * math.sin(2*alpha)
    CA_conicfrostrum = (1-(NoseRadius/BaseRadius)**2 * (math.cos(delta)**2)) * (2 * (math.sin(delta)**2) * (math.cos(alpha)**2) + (math.cos(delta)**2)*(math.sin(delta)**2))

    CA_body = (NoseRadius/BaseRadius)**2 * CD_nose + CA_conicfrostrum
    CN_body = (NoseRadius/BaseRadius)**2 * CL_nose + CN_conicfrostrum

    CL_body = CN_body * math.cos(alpha) - CA_body * math.sin(alpha)
    CD_body = CA_body * math.cos(alpha) + CN_body * math.sin(alpha)
    return CL_body, CD_body


