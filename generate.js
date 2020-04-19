fs = require('fs');

var Integrators = {
  Euler    : [[1]],
  Midpoint : [[.5,.5],[0, 1]],
  Heun     : [[1, 1],[.5,.5]],
  Ralston  : [[2/3,2/3],[.25,.75]],
  K3       : [[.5,.5],[1,-1,2],[1/6,2/3,1/6]],
  SSP33    : [[1,1],[.5,.25,.25],[1/6,1/6,2/3]],
  SSP43    : [[.5,.5],[1,.5,.5],[.5,1/6,1/6,1/6],[1/6,1/6,1/6,1/2]],
  RK4      : [[.5,.5],[.5,0,.5],[1,0,0,1],[1/6,1/3,1/3,1/6]],
  RK38     : [[1/3,1/3],[2/3,-1/3,1],[1,1,-1,1],[1/8,3/8,3/8,1/8]]
};


// f is a func of time t and state y
// y is the initial state, t is the time, h is the timestep
// updated y is returned.
var integrate=(m,f,y,t,h)=>{
  for (var k=[],ki=0; ki<m.length; ki++) {
    var _y=y.slice(), dt=ki?((m[ki-1][0])*h):0;
    for (var l=0; l<_y.length; l++) for (var j=1; j<=ki; j++) _y[l]=_y[l]+h*(m[ki-1][j])*(k[ki-1][l]);
    k[ki]=f(t+dt,_y,dt); 
  }
  for (var r=y.slice(),l=0; l<_y.length; l++) for (var j=0; j<k.length; j++) r[l]=r[l]+h*(k[j][l])*(m[ki-1][j]);
  return r;
}


function get_solution(dt, N, I0, R0, D_incbation, D_infectious, D_recovery_mild, D_hospital_lag, D_recovery_severe, D_death, P_SEVERE, CFR, InterventionTime, InterventionAmt, duration) {
  var interpolation_steps = 40
  var steps = 110*interpolation_steps
  var dt = dt/interpolation_steps
  var sample_step = interpolation_steps
  var method = Integrators["RK4"]
  function f(t, x){
    // SEIR ODE
    if (t > InterventionTime && t < InterventionTime + duration){
      var beta = (InterventionAmt)*R0/(D_infectious)
    } else if (t > InterventionTime + duration) {
      var beta = 0.5*R0/(D_infectious)        
    } else {
      var beta = R0/(D_infectious)
    }
    var a     = 1/D_incbation
    var gamma = 1/D_infectious
    
    var S        = x[0] // Susectable
    var E        = x[1] // Exposed
    var I        = x[2] // Infectious 
    var Mild     = x[3] // Recovering (Mild)     
    var Severe   = x[4] // Recovering (Severe at home)
    var Severe_H = x[5] // Recovering (Severe in hospital)
    var Fatal    = x[6] // Recovering (Fatal)
    var R_Mild   = x[7] // Recovered
    var R_Severe = x[8] // Recovered
    var R_Fatal  = x[9] // Dead
    var p_severe = P_SEVERE
    var p_fatal  = CFR
    var p_mild   = 1 - P_SEVERE - CFR
    var dS        = -beta*I*S
    var dE        =  beta*I*S - a*E
    var dI        =  a*E - gamma*I
    var dMild     =  p_mild*gamma*I   - (1/D_recovery_mild)*Mild
    var dSevere   =  p_severe*gamma*I - (1/D_hospital_lag)*Severe
    var dSevere_H =  (1/D_hospital_lag)*Severe - (1/D_recovery_severe)*Severe_H
    var dFatal    =  p_fatal*gamma*I  - (1/D_death)*Fatal
    var dR_Mild   =  (1/D_recovery_mild)*Mild
    var dR_Severe =  (1/D_recovery_severe)*Severe_H
    var dR_Fatal  =  (1/D_death)*Fatal
    //      0   1   2   3      4        5          6       7        8          9
    return [dS, dE, dI, dMild, dSevere, dSevere_H, dFatal, dR_Mild, dR_Severe, dR_Fatal]
  }
  var v = [1 - I0/N, 0, I0/N, 0, 0, 0, 0, 0, 0, 0]
  var t = 0
  var P  = []
  var TI = []
  var Iters = []
  while (steps--) { 
    if ((steps+1) % (sample_step) == 0) {
          //    Dead   Hospital          Recovered        Infectious   Exposed
      P.push([ N*v[9], N*(v[5]+v[6]),  N*(v[7] + v[8]), N*v[2],    N*v[1] ])
      Iters.push(v)
      TI.push(N*(1-v[0]))
      // console.log((v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7] + v[8] + v[9]))
      // console.log(v[0] , v[1] , v[2] , v[3] , v[4] , v[5] , v[6] , v[7] , v[8] , v[9])
    }
    v =integrate(method,f,v,t,dt); 
    t+=dt
  }
  return {"P": P, 
          "deaths": N*v[6], 
          "total": 1-v[0],
          "total_infected": TI,
          "Iters":Iters,
          "dIters": f}
}


let Time_to_death     = 32;

let I0                = 1;
let R0                = 2.2;
let D_incbation       = 5.2;
let D_infectious      = 2.9;
let D_recovery_mild   = (14 - 2.9);
let D_recovery_severe = (31.5 - 2.9);
let D_hospital_lag    = 5;
let D_death           = Time_to_death - D_infectious;
let CFR               = 0.02;
let InterventionTime  = 100;
let OMInterventionAmt = 2/3;
let InterventionAmt   = 1 - OMInterventionAmt;
let Time              = 220;
let Xmax              = 110000;
let dt                = 2;
let P_SEVERE          = 0.2;
let duration          = 7*12*1e10;

var logger = fs.createWriteStream('out.csv');


for (let population = 20000.; population < 10e6; population*=1.07) {
    for (let InterventionTime = 20; InterventionTime < 200; InterventionTime += 5) {

        // let population = 7e6;
        // let InterventionTime  = 100;

        let logN              = Math.log(population);
        let N                 = Math.exp(logN);

        let line = "";
        
        line += population;
        line +=  "," + InterventionTime;

        r = get_solution(dt, N, I0, R0, D_incbation, D_infectious,
                         D_recovery_mild, D_hospital_lag, D_recovery_severe, 
                         D_death, P_SEVERE, CFR, InterventionTime, 
                         InterventionAmt, duration);

        r["P"].forEach((x) => {
            line += "," + x[0];  // Dead
        });

        r["P"].forEach((x) => {
            line += "," + x[0];  // Hospital
        });

        logger.write(line + "\n");

}
}

logger.end();
