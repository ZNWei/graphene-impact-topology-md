"""Unmodified logistic_fit function extracted from archived round-2 analysis.
The only changes are the module imports and this provenance note.
"""
from __future__ import annotations
import math
from typing import Any

def logistic_fit(summary: list[dict[str, Any]], cfg: dict[str, Any]) -> dict[str, Any]:
    """Finite ordinary grouped-binomial logistic fit for three or more speed groups."""
    model = cfg["analysis_model"]; center=float(model["velocity_center_mps"]); scale=float(model["velocity_scale_mps"])
    if len(summary) < 3 or any(int(x["n"]) != 10 for x in summary): return {"status":"not_estimable_incomplete_design","reason":"at least three complete ten-replicate velocity groups are required"}
    probs=[float(x["event_probability"]) for x in summary]
    if any(b<a for a,b in zip(probs,probs[1:])): return {"status":"not_estimable_nonmonotonic","reason":"observed event fractions decrease with velocity"}
    if not (probs[0]<0.5<probs[-1]): return {"status":"not_estimable_unbracketed","reason":"0.5 event probability is not bracketed"}
    mixed=[x for x in summary if 0<int(x["events"])<int(x["n"])]
    if len(mixed)<2: return {"status":"not_estimable_monotone_likelihood","reason":"fewer than two distinct non-extreme velocity groups leave the ordinary grouped-logistic likelihood unbounded"}
    x=[(float(v["velocity_mps"])-center)/scale for v in summary]; y=[int(v["events"]) for v in summary]; n=[int(v["n"]) for v in summary]
    b0,b1=0.0,1.0
    for _ in range(int(model["maximum_iterations"])):
        eta=[b0+b1*xi for xi in x]; p=[1/(1+math.exp(-max(-35,min(35,e)))) for e in eta]; w=[ni*pi*(1-pi) for ni,pi in zip(n,p)]
        g0=sum(yi-ni*pi for yi,ni,pi in zip(y,n,p)); g1=sum(xi*(yi-ni*pi) for xi,yi,ni,pi in zip(x,y,n,p)); h00=sum(w); h01=sum(wi*xi for wi,xi in zip(w,x)); h11=sum(wi*xi*xi for wi,xi in zip(w,x)); det=h00*h11-h01*h01
        if det<=1e-14: return {"status":"not_estimable_singular_information","reason":"logistic information matrix is singular"}
        d0=(h11*g0-h01*g1)/det; d1=(-h01*g0+h00*g1)/det; b0+=d0; b1+=d1
        if max(abs(d0),abs(d1))<float(model["convergence_tolerance"]): break
    else: return {"status":"not_estimable_nonconvergent","reason":"finite grouped-logistic fit did not converge"}
    if b1<=0: return {"status":"not_estimable_nonpositive_slope","reason":"fitted slope is not positive"}
    eta=[b0+b1*xi for xi in x]; p=[1/(1+math.exp(-max(-35,min(35,e)))) for e in eta]; w=[ni*pi*(1-pi) for ni,pi in zip(n,p)]; h00=sum(w); h01=sum(wi*xi for wi,xi in zip(w,x)); h11=sum(wi*xi*xi for wi,xi in zip(w,x)); det=h00*h11-h01*h01
    if det<=1e-14: return {"status":"not_estimable_singular_information","reason":"final information matrix is singular"}
    vt=center-scale*b0/b1; low=min(x["velocity_mps"] for x in summary); high=max(x["velocity_mps"] for x in summary)
    if not low<=vt<=high: return {"status":"not_estimable_threshold_outside_range","reason":"fitted V_T,50 lies outside tested range"}
    cov00,cov01,cov11=h11/det,-h01/det,h00/det; g0=-scale/b1; g1=scale*b0/(b1*b1); se=math.sqrt(max(0,g0*g0*cov00+2*g0*g1*cov01+g1*g1*cov11))
    return {"status":"estimated","estimator":"pooled grouped binomial logistic regression","slope_per_1000_mps":b1,"v_t_50_mps":vt,"wald_95_ci_mps":[vt-1.959963984540054*se,vt+1.959963984540054*se],"warning":"This is a fragmentation-mediated projectile-origin mass-transmission threshold, not conventional ballistic V50."}
