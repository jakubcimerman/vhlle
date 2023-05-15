class EoS;
class EoSaux;
class EoSChiral : public EoS {
private:
 EoSaux *eossmall, *eosbig;

public:
 EoSChiral(void);
 ~EoSChiral(void);

 // eos_full returns EoS interpolated from the small EoS table, or 
 // from a big EoS table, or extrapolated EoS otherwise
 void eos_full(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 // final EoS function: calls eos_full() and re-evaluates EoS at lower
 // baryon density if the evaluated pressure is too small
 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 // p_full returns pressure interpolated from the small EoS table, or 
 // from a big EoS table, or extrapolated pressure otherwise
 double p_full(double e, double nb, double nq, double ns);
 // final EoS function: calls p_full() and re-evaluates pressure at lower
 // baryon density if the evaluated pressure is too small
 virtual double p(double e, double nb, double nq, double ns);
};
