MODULE phys
IMPLICIT NONE
integer, parameter :: dp = kind(1.0d0)
real(kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
!-------------Basic physical constants-------------------
!
!The following five are the most accurate 1986 values
REAL (KIND=8), PARAMETER :: h = 6.626075540D-34    !Planck's constant
REAL (KIND=8), PARAMETER :: c = 2.99792458D+08       !Speed of light in a vacuum
REAL (KIND=8), PARAMETER :: c2 = c**2
REAL (KIND=8), PARAMETER :: k = 1.38065812D-23     !Boltzman thermodynamic constant
REAL (KIND=8), PARAMETER :: kb = 1.36D-22          !Boltzmann constant in atm.molec^-1.cm^3.K^-1
REAL (KIND=8), PARAMETER :: sigma = 5.67051196D-08  !Stefan-Boltzman constant
REAL (KIND=8), PARAMETER :: G_cst = 6.67428D-11        !Gravitational constant (2006 measurements)

!-----------Thermodynamic constants------------
REAL (KIND=8), PARAMETER :: N_avogadro = 6.022136736D+23  !Avogadro's number
REAL (KIND=8), PARAMETER :: Rstar = 1000D0*k*N_avogadro   !Universal gas constant

!Set properties of individual gases

!Thermodynamic properties of dry Earth air
REAL (KIND=8), PARAMETER :: cp_dry_Earth_air = 1004D0
REAL (KIND=8), PARAMETER :: MolecularWeight_dry_Earth_air = 28.97D0
REAL (KIND=8), PARAMETER :: gamma_dry_Earth_air = 1.4003D0
!--------------------------------------------------------    
!H2O
REAL (KIND=8), PARAMETER :: H2O_CriticalPointT = 6.471000D+02
REAL (KIND=8), PARAMETER :: H2O_CriticalPointP = 2.210000D+07
REAL (KIND=8), PARAMETER :: H2O_TriplePointT = 2.731500D+02
REAL (KIND=8), PARAMETER :: H2O_TriplePointP = 6.110000D+02
REAL (KIND=8), PARAMETER :: H2O_L_vaporization_BoilingPoint = 2.255000D+06
REAL (KIND=8), PARAMETER :: H2O_L_vaporization_TriplePoint = 2.493000D+06
REAL (KIND=8), PARAMETER :: H2O_L_fusion = 3.340000D+05
REAL (KIND=8), PARAMETER :: H2O_L_sublimation = 2.840000D+06
REAL (KIND=8), PARAMETER :: H2O_rho_liquid_BoilingPoint = 9.584000D+02
REAL (KIND=8), PARAMETER :: H2O_rho_liquid_TriplePoint = 9.998700D+02
REAL (KIND=8), PARAMETER :: H2O_rho_solid = 9.170000D+02
REAL (KIND=8), PARAMETER :: H2O_cp = 1.847000D+03
REAL (KIND=8), PARAMETER :: H2O_gamma = 1.331000D+00
REAL (KIND=8), PARAMETER :: H2O_MolecularWeight = 1.800000D+01
CHARACTER (LEN = 5), PARAMETER :: H2O_name = 'Water'
CHARACTER (LEN = 3), PARAMETER :: H2O_formula = 'H2O'
REAL (KIND=8), PARAMETER :: H2O_L_vaporization=2.493000D+06
REAL (KIND=8), PARAMETER :: H2O_rho_liquid=9.998700D+02
!--------------------------------------------------------
!CH4
REAL (KIND=8), PARAMETER :: CH4_CriticalPointT = 1.904400D+02
REAL (KIND=8), PARAMETER :: CH4_CriticalPointP = 4.596000D+06
REAL (KIND=8), PARAMETER :: CH4_TriplePointT = 9.067000D+01
REAL (KIND=8), PARAMETER :: CH4_TriplePointP = 1.170000D+04
REAL (KIND=8), PARAMETER :: CH4_L_vaporization_BoilingPoint = 5.100000D+05
REAL (KIND=8), PARAMETER :: CH4_L_vaporization_TriplePoint = 5.360000D+05
REAL (KIND=8), PARAMETER :: CH4_L_fusion = 5.868000D+04
REAL (KIND=8), PARAMETER :: CH4_L_sublimation = 5.950000D+05
REAL (KIND=8), PARAMETER :: CH4_rho_liquid_BoilingPoint = 4.502000D+02
REAL (KIND=8), PARAMETER :: CH4_rho_liquid_TriplePoint = 0D0
REAL (KIND=8), PARAMETER :: CH4_rho_solid = 5.093000D+02
REAL (KIND=8), PARAMETER :: CH4_cp = 2.195000D+03
REAL (KIND=8), PARAMETER :: CH4_gamma = 1.305000D+00
REAL (KIND=8), PARAMETER :: CH4_MolecularWeight = 1.600000D+01
CHARACTER (LEN = 7), PARAMETER :: CH4_name = 'Methane'
CHARACTER (LEN = 3), PARAMETER :: CH4_formula = 'CH4'
REAL (KIND=8), PARAMETER :: CH4_L_vaporization=5.360000D+05
REAL (KIND=8), PARAMETER :: CH4_rho_liquid=4.502000D+02
!--------------------------------------------------------
!CO2
REAL (KIND=8), PARAMETER :: CO2_CriticalPointT = 3.042000D+02
REAL (KIND=8), PARAMETER :: CO2_CriticalPointP = 7.382500D+06
REAL (KIND=8), PARAMETER :: CO2_TriplePointT = 2.165400D+02
REAL (KIND=8), PARAMETER :: CO2_TriplePointP = 5.185000D+05
REAL (KIND=8), PARAMETER :: CO2_L_vaporization_BoilingPoint = 0D0
REAL (KIND=8), PARAMETER :: CO2_L_vaporization_TriplePoint = 3.970000D+05
REAL (KIND=8), PARAMETER :: CO2_L_fusion = 1.960000D+05
REAL (KIND=8), PARAMETER :: CO2_L_sublimation = 5.930000D+05
REAL (KIND=8), PARAMETER :: CO2_rho_liquid_BoilingPoint = 1.032000D+03
REAL (KIND=8), PARAMETER :: CO2_rho_liquid_TriplePoint = 1.110000D+03
REAL (KIND=8), PARAMETER :: CO2_rho_solid = 1.562000D+03
REAL (KIND=8), PARAMETER :: CO2_cp = 8.200000D+02
REAL (KIND=8), PARAMETER :: CO2_gamma = 1.294000D+00
REAL (KIND=8), PARAMETER :: CO2_MolecularWeight = 4.400000D+01
CHARACTER (LEN = 14), PARAMETER :: CO2_name = 'Carbon Dioxide'
CHARACTER (LEN = 3), PARAMETER :: CO2_formula = 'CO2'
REAL (KIND=8), PARAMETER :: CO2_L_vaporization=3.970000D+05
REAL (KIND=8), PARAMETER :: CO2_rho_liquid=1.110000D+03
!--------------------------------------------------------
!CO
REAL (KIND=8), PARAMETER :: CO_CriticalPointT = 1.134450D+02
REAL (KIND=8), PARAMETER :: CO_CriticalPointP = 3.498750D+06
REAL (KIND=8), PARAMETER :: CO_TriplePointT = 6.795000D+01
REAL (KIND=8), PARAMETER :: CO_TriplePointP = 1.530000D+04
REAL (KIND=8), PARAMETER :: CO_L_vaporization_BoilingPoint = 0D0
REAL (KIND=8), PARAMETER :: CO_L_vaporization_TriplePoint = 2.142857D+05
REAL (KIND=8), PARAMETER :: CO_L_fusion = 0D0
REAL (KIND=8), PARAMETER :: CO_L_sublimation = 2.7142857D+05
REAL (KIND=8), PARAMETER :: CO_rho_liquid_BoilingPoint = 0D0
REAL (KIND=8), PARAMETER :: CO_rho_liquid_TriplePoint = 0D0
REAL (KIND=8), PARAMETER :: CO_rho_solid = 0D0
REAL (KIND=8), PARAMETER :: CO_cp = 1.04000D+03
REAL (KIND=8), PARAMETER :: CO_gamma = 50.241546D+00
REAL (KIND=8), PARAMETER :: CO_MolecularWeight = 2.800970D+01
CHARACTER (LEN = 14), PARAMETER :: CO_name = 'Carbon Monoxide'
CHARACTER (LEN = 3), PARAMETER :: CO_formula = 'CO'
REAL (KIND=8), PARAMETER :: CO_L_vaporization=2.142857D+05
REAL (KIND=8), PARAMETER :: CO_rho_liquid = 0D0
!--------------------------------------------------------
!N2
REAL (KIND=8), PARAMETER :: N2_CriticalPointT = 1.262000D+02
REAL (KIND=8), PARAMETER :: N2_CriticalPointP = 3.400000D+06
REAL (KIND=8), PARAMETER :: N2_TriplePointT = 6.314000D+01
REAL (KIND=8), PARAMETER :: N2_TriplePointP = 1.253000D+04
REAL (KIND=8), PARAMETER :: N2_L_vaporization_BoilingPoint = 1.980000D+05
REAL (KIND=8), PARAMETER :: N2_L_vaporization_TriplePoint = 2.180000D+05
REAL (KIND=8), PARAMETER :: N2_L_fusion = 2.573000D+04
REAL (KIND=8), PARAMETER :: N2_L_sublimation = 2.437000D+05
REAL (KIND=8), PARAMETER :: N2_rho_liquid_BoilingPoint = 8.086000D+02
REAL (KIND=8), PARAMETER :: N2_rho_liquid_TriplePoint = 0D0
REAL (KIND=8), PARAMETER :: N2_rho_solid = 1.026000D+03
REAL (KIND=8), PARAMETER :: N2_cp = 1.037000D+03
REAL (KIND=8), PARAMETER :: N2_gamma = 1.403000D+00
REAL (KIND=8), PARAMETER :: N2_MolecularWeight = 2.800000D+01
CHARACTER (LEN = 8), PARAMETER :: N2_name = 'Nitrogen'
CHARACTER (LEN = 2), PARAMETER :: N2_formula = 'N2'
REAL (KIND=8), PARAMETER :: N2_L_vaporization=2.180000D+05
REAL (KIND=8), PARAMETER :: N2_rho_liquid=8.086000D+02
!--------------------------------------------------------
!O2
REAL (KIND=8), PARAMETER :: O2_CriticalPointT = 1.545400D+02
REAL (KIND=8), PARAMETER :: O2_CriticalPointP = 5.043000D+06
REAL (KIND=8), PARAMETER :: O2_TriplePointT = 5.430000D+01
REAL (KIND=8), PARAMETER :: O2_TriplePointP = 1.500000D+02
REAL (KIND=8), PARAMETER :: O2_L_vaporization_BoilingPoint = 2.130000D+05
REAL (KIND=8), PARAMETER :: O2_L_vaporization_TriplePoint = 2.420000D+05
REAL (KIND=8), PARAMETER :: O2_L_fusion = 1.390000D+04
REAL (KIND=8), PARAMETER :: O2_L_sublimation = 2.560000D+05
REAL (KIND=8), PARAMETER :: O2_rho_liquid_BoilingPoint = 1.141000D+03
REAL (KIND=8), PARAMETER :: O2_rho_liquid_TriplePoint = 1.307000D+03
REAL (KIND=8), PARAMETER :: O2_rho_solid = 1.351000D+03
REAL (KIND=8), PARAMETER :: O2_cp = 9.160000D+02
REAL (KIND=8), PARAMETER :: O2_gamma = 1.393000D+00
REAL (KIND=8), PARAMETER :: O2_MolecularWeight = 3.200000D+01
CHARACTER (LEN = 6), PARAMETER :: O2_name = 'Oxygen'
CHARACTER (LEN = 2), PARAMETER :: O2_formula = 'O2'
REAL (KIND=8), PARAMETER :: O2_L_vaporization=2.420000D+05
REAL (KIND=8), PARAMETER :: O2_rho_liquid=1.307000D+03
!--------------------------------------------------------
!H2
REAL (KIND=8), PARAMETER :: H2_CriticalPointT = 3.320000D+01
REAL (KIND=8), PARAMETER :: H2_CriticalPointP = 1.298000D+06
REAL (KIND=8), PARAMETER :: H2_TriplePointT = 1.395000D+01
REAL (KIND=8), PARAMETER :: H2_TriplePointP = 7.200000D+03
REAL (KIND=8), PARAMETER :: H2_L_vaporization_BoilingPoint = 4.540000D+05
REAL (KIND=8), PARAMETER :: H2_L_vaporization_TriplePoint = 0D0
REAL (KIND=8), PARAMETER :: H2_L_fusion = 5.820000D+04
REAL (KIND=8), PARAMETER :: H2_L_sublimation = 0D0
REAL (KIND=8), PARAMETER :: H2_rho_liquid_BoilingPoint = 7.097000D+01
REAL (KIND=8), PARAMETER :: H2_rho_liquid_TriplePoint = 0D0
REAL (KIND=8), PARAMETER :: H2_rho_solid = 8.800000D+01
REAL (KIND=8), PARAMETER :: H2_cp = 1.423000D+04
REAL (KIND=8), PARAMETER :: H2_gamma = 1.384000D+00
REAL (KIND=8), PARAMETER :: H2_MolecularWeight = 2.000000D+00
CHARACTER (LEN = 8), PARAMETER :: H2_name = 'Hydrogen'
CHARACTER (LEN = 2), PARAMETER :: H2_formula = 'H2'
REAL (KIND=8), PARAMETER :: H2_L_vaporization=4.540000D+05
REAL (KIND=8), PARAMETER :: H2_rho_liquid=7.097000D+01
!--------------------------------------------------------
!He
REAL (KIND=8), PARAMETER :: He_CriticalPointT = 5.100000D+00
REAL (KIND=8), PARAMETER :: He_CriticalPointP = 2.280000D+05
REAL (KIND=8), PARAMETER :: He_TriplePointT = 2.170000D+00
REAL (KIND=8), PARAMETER :: He_TriplePointP = 5.070000D+03
REAL (KIND=8), PARAMETER :: He_L_vaporization_BoilingPoint = 2.030000D+04
REAL (KIND=8), PARAMETER :: He_L_vaporization_TriplePoint = 0D0
REAL (KIND=8), PARAMETER :: He_L_fusion = 0D0
REAL (KIND=8), PARAMETER :: He_L_sublimation = 0D0
REAL (KIND=8), PARAMETER :: He_rho_liquid_BoilingPoint = 1.249600D0+02
REAL (KIND=8), PARAMETER :: He_rho_liquid_TriplePoint = 0D0
REAL (KIND=8), PARAMETER :: He_rho_solid = 2.000000D+02
REAL (KIND=8), PARAMETER :: He_cp = 5.196000D+03
REAL (KIND=8), PARAMETER :: He_gamma = 1.664000D+00
REAL (KIND=8), PARAMETER :: He_MolecularWeight = 4.000000D+00
CHARACTER (LEN = 6), PARAMETER :: He_name = 'Helium'
CHARACTER (LEN = 2), PARAMETER :: He_formula = 'He'
REAL (KIND=8), PARAMETER :: He_L_vaporization=2.030000D+04
REAL (KIND=8), PARAMETER :: He_rho_liquid=1.249600D+02
!--------------------------------------------------------
!NH3
REAL (KIND=8), PARAMETER :: NH3_CriticalPointT = 4.055000D+02
REAL (KIND=8), PARAMETER :: NH3_CriticalPointP = 1.128000D+07
REAL (KIND=8), PARAMETER :: NH3_TriplePointT = 1.954000D+02
REAL (KIND=8), PARAMETER :: NH3_TriplePointP = 6.100000D+03
REAL (KIND=8), PARAMETER :: NH3_L_vaporization_BoilingPoint = 1.371000D+06
REAL (KIND=8), PARAMETER :: NH3_L_vaporization_TriplePoint = 1.658000D+06
REAL (KIND=8), PARAMETER :: NH3_L_fusion = 3.314000D+05
REAL (KIND=8), PARAMETER :: NH3_L_sublimation = 1.989000D+06
REAL (KIND=8), PARAMETER :: NH3_rho_liquid_BoilingPoint = 6.820000D+02
REAL (KIND=8), PARAMETER :: NH3_rho_liquid_TriplePoint = 7.342000D+02
REAL (KIND=8), PARAMETER :: NH3_rho_solid = 8.226000D+02
REAL (KIND=8), PARAMETER :: NH3_cp = 2.060000D+03
REAL (KIND=8), PARAMETER :: NH3_gamma = 1.309000D+00
REAL (KIND=8), PARAMETER :: NH3_MolecularWeight = 1.700000D+01
CHARACTER (LEN = 7), PARAMETER :: NH3_name = 'Ammonia'
CHARACTER (LEN = 3), PARAMETER :: NH3_formula = 'NH3'
REAL (KIND=8), PARAMETER :: NH3_L_vaporization=1.658000D+06
REAL (KIND=8), PARAMETER :: NH3_rho_liquid=7.342000D+02

CONTAINS
!Planck function (of frequency)
    FUNCTION B(nu,T) RESULT(res)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN) :: nu,T
        REAL (KIND=8) :: u,res
        u=min(h*nu/(k*T),200D0)
        res=(2D0*h*nu**3/c**2)/(exp(u)-1D0)
    END FUNCTION B
    
    FUNCTION Blam(lam,T) RESULT(res)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN) :: lam,T
        REAL (KIND=8) :: u,res
        u=h*c/(lam*k*T)
        res=(2D0*h*c**2/lam**5)/(exp(u)-1D0)
    END FUNCTION Blam

    FUNCTION dB(nu,T) RESULT(res)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN) :: nu,T
        REAL (KIND=8) :: u,res
        u =min(h*nu/(k*T),200D0)
        res=(2D0*h**2*nu**4/(k*c**2*T**2))*(exp(u)/(exp(u)-1D0)**2)
        !print*, T
    END FUNCTION dB
    
    !Saturation vapor pressure over ice (Smithsonian formula)
    !Input: Kelvin. Output: Pascal
    FUNCTION satvpi(T) RESULT(res)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN) :: T
        REAL (KIND=8) :: esbasi,tbasi,aa,b,c,e,esice,res
        !Compute es over ice (valid between -153 c and 0 c)
        !see smithsonian meteorological tables page 350
   
        !Original source: GFDL climate model, circa 1995
        esbasi = 6107.1D0
        tbasi  = 273.16D0
   
        aa     = -9.09718D0 *(tbasi/T-1.0D0)
        b      = -3.56654D0 *log10(tbasi/T)
        c      = 0.876793D0*(1.0D0-T/tbasi)
        e      = log10(esbasi)
        esice  = 10D0**(aa+b+c+e)
        res =.1D0*esice  !Convert to Pascals
    END FUNCTION satvpi

    !Saturation vapor pressure over liquid water (Smithsonian formula)
    !Input: Kelvin. Output: Pascal
    FUNCTION satvpw(T) RESULT(res)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN) :: T
        REAL (KIND=8) :: esbasw,tbasw,aa,b,c,d,e,esh2O,res
        !compute es over liquid water between -20c and freezing.
        !see smithsonian meteorological tables page 350.
   
        !Original source: GFDL climate model, circa 1995
        esbasw = 1013246.0D0
        tbasw  = 373.16D0

        aa     = -7.90298D0*(tbasw/T-1D0)
        b      =  5.02808D0*log10(tbasw/T)
        c      = -1.3816D-07*(  10D0**( ((1D0-T/tbasw)*11.344D0)-1D0 )  )
        d      =  8.1328D-03*(  10D0**( ((tbasw/T-1D0)*(-3.49149D0))-1D0)  )
        e      = log10(esbasw)
        esh2O  = 10D0**(aa+b+c+d+e)
        res = .1D0*esh2O  !Convert to Pascals
    END FUNCTION satvpw

    ! An alternate formula for saturation vapor pressure over liquid water
    FUNCTION satvpw_Heymsfield(T) RESULT(res)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN) :: T
        REAL (KIND=8) :: ts,sr,ar,br,cr,dw,er,vp,res
        ts = 373.16D0
        sr = 3.0057166D0
        ! Vapor pressure over water. Heymsfield formula
        ar = ts/T
        br = 7.90298D0*(ar-1D0)
        cr  = 5.02808D0*log10(ar);
        dw = (1.3816D-07)*(10D0**(11.344D0*(1D0-1D0/ar))-1D0)
        er = 8.1328D-03*((10D0**(-(3.49149D0*(ar-1D0))) )-1D0)
        vp = 10D0**(cr-dw+er+sr-br)
        vp = vp*1.0D+02
        res = vp
    END FUNCTION satvpw_Heymsfield
      
    FUNCTION  satvpg(T) RESULT(res)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN) :: T
        REAL (KIND=8) :: res
         
        !This is the saturation vapor pressure computation used in the
        !GFDL climate model.  It blends over from water saturation to
        !ice saturation as the temperature falls below 0C.
        IF ((T-273.16D0) .LT. -20D0) THEN
            res=satvpi(T)
        ENDIF
        IF ( ((T-273.16D0) .GE. -20D0) .AND. ((T-273.16D0) .LE. 0D0)) THEN
            res=0.05*(273.16D0-T)*satvpi(T) + 0.05D0*(T-253.16D0)*satvpw(T)
        ENDIF
        IF ((T-273.16D0) .GT. 0D0) THEN
            res=satvpw(T)
        ENDIF
    END FUNCTION satvpg

    !Saturation vapor pressure for any substance, computed using
    !the simplified form of Clausius-Clapeyron assuming the perfect
    !gas law and constant latent heat
    FUNCTION satvps(T,T0,e0,MolecularWeight,LatentHeat) RESULT(res)
        IMPLICIT NONE
        REAL (KIND=8), INTENT(IN) :: T,T0,e0,MolecularWeight,LatentHeat
        REAL (KIND=8) :: Rv,res
        Rv = Rstar/MolecularWeight 
        res = e0*exp(-(LatentHeat/Rv)*(1D0/T - 1D0/T0))
    END FUNCTION satvps
END MODULE phys
