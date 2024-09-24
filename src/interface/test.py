import pymc
import cronos
import numpy as np

def ode_define():

  # Define DAG
  DAG = pymc.FFGraph()
  X1  = pymc.FFVar(DAG,"X1")
  X2  = pymc.FFVar(DAG,"X2")
  P   = pymc.FFVar(DAG,"P" )
  X10 = pymc.FFVar(DAG,"X10" )
#  X10.set(1.2)

  # Define IVP
  ODE = cronos.ODESLV()
  ODE.set_dag( DAG )
  ODE.set_time( [ 0., 5., 10. ] )
#  print( ODE.val_stage, ODE.var_time )
  ODE.set_state( [X1,X2] )
#  print( ODE.var_state )
  ODE.set_parameter( [P] )
#  print( ODE.var_parameter )
  ODE.set_constant( [X10] )
#  print( ODE.var_constant )
  ODE.set_differential( [ P*X1*(1-X2), P*X2*(X1-1) ] )
#  print( ODE.eqn_differential )
  ODE.set_initial( [ X10, 1.1+0.01*(P-3) ] )
#  print( ODE.eqn_initial )
  ODE.set_function( [ X1*X2, P*pymc.sqr(X1) ] )
#  print( ODE.eqn_function )

  ODE.options.DISPLEVEL = 1
  ODE.options.INTMETH   = ODE.options.MSBDF
  ODE.options.NLINSOL   = ODE.options.NEWTON #FIXEDPOINT
  ODE.options.LINSOL    = ODE.options.DENSE  #DIAG
  ODE.setup()

  stat = ODE.solve_state( [2.95], [1.2] )
  stat = ODE.solve_state( [2.95], [1.1] )

  print( "status:", stat )
  print( "final time:", ODE.final_time )
  print( "final stage:", ODE.final_stage )
  print( ODE.val_state )

  ODE.solve_sensitivity( [2.95], [1.1] )
  print( ODE.val_function_gradient )

  ODE.solve_adjoint( [2.95], [1.1] )
  print( ODE.val_function_gradient )

  return ODE

  
def ode_copy( ODE ):

  ODE_copy = cronos.ODESLV()
  ODE_copy.set( ODE )
  ODE_copy.options = ODE.options
  ODE_copy.options.RESRECORD = 100
  ODE_copy.setup()

  stat = ODE_copy.solve_sensitivity( [2.95], [1.1] )
  print( ODE_copy.results_state, ODE_copy.results_sensitivity )

  return ODE_copy

  
def ode_diff( ODE ):

  ODE_diff = ODE.fdiff( ODE.var_parameter+ODE.var_constant ) #parameter )
  ODE_diff.setup()
  
  stat = ODE.solve_state( [2.95], [1.1] )
  stat = ODE_diff.solve_state( [2.95], [1.1] )

  return ODE_diff

  
def ode_dag( ODE ):

  OpODE = cronos.FFODE()
  DAG = pymc.FFGraph()
  P = pymc.FFVar( DAG, "P" )
  C = pymc.FFVar( DAG, "C" )
  F = OpODE( [P], [C], ODE )
  SGF = DAG.subgraph( F )
  DAG.output( SGF )
  DAG.dot_script( F, "F.dot" )

  print( "F @(2.95,1.1): ", DAG.eval( F, [P,C], [2.95,1.1] ) )

  DFDP = DAG.bdiff( F, [P] )
  SGF = DAG.subgraph( DFDP[2] )
  DAG.output( SGF )
  DAG.dot_script( DFDP[2], "DFDP.dot" )

  print( "DFDP @(2.95,1.1): ", DAG.eval( DFDP[2], [P,C], [2.95,1.1] ) )
  
  OpODE.options.DIFF = OpODE.options.SYM_PC
  DFDP = DAG.bdiff( F, [P,C] )
  SGF = DAG.subgraph( DFDP[2] )
  DAG.output( SGF )
  DAG.dot_script( DFDP[2], "DFDP.dot" )

  print( "DFDP @(2.95,1.1): ", DAG.eval( DFDP[2], [P,C], [2.95,1.1] ) )
  
def ode_dag2( ODE ):

  OpODE = cronos.FFODE()
  DAG = pymc.FFGraph()
  P1 = pymc.FFVar( DAG, "P1" )
  P2 = pymc.FFVar( DAG, "P2" )
  C  = pymc.FFVar( DAG, "C" )
  F1 = OpODE( [P1], [C], ODE )
  F2 = OpODE( [P2], [C], ODE )
  F12 = F1+F2
  SGF = DAG.subgraph( F12 )
  DAG.output( SGF )
  DAG.dot_script( F12, "F12.dot" )

  print( "F12 @(2.95,3.05,1.1): ", DAG.eval( F12, [P1,P2,C], [2.95,3.05,1.1] ) )

  DF12DP = DAG.bdiff( F12, [P1,P2] )
  SGF12 = DAG.subgraph( DF12DP[2] )
  DAG.output( SGF12 )
  DAG.dot_script( DF12DP[2], "DF12DP.dot" )

  print( "DF12DP @(2.95,3.05,1.1): ", DAG.eval( DF12DP[2], [P1,P2,C], [2.95,3.05,1.1] ) )


def ode2_define( NS ):

  # Define DAG
  DAG = pymc.FFGraph()
  X = [pymc.FFVar(DAG,"X")]
  Q = [pymc.FFVar(DAG,"Q")]
  U    = []
  RHS  = []
  QUAD = []
  F    = []
  for i in range(NS):
    U.append( pymc.FFVar(DAG,"U"+str(i)) )
    RHS.append( [ U[i] - X[0] ] )
    QUAD.append( [ 0.5 * pymc.sqr( U[i] ) ] )
    F.append( [ Q[0], pymc.FFVar(0.) ] )
  F[NS-1][1] = X[0]
   
  # Define IVP
  ODE = cronos.ODESLV()
  ODE.set_dag( DAG )
  ODE.set_time( np.arange(0., 1.01, 1./NS).tolist() )
  print( ODE.val_stage, ODE.var_time )
  ODE.set_state( X )
  ODE.set_parameter( U )
  ODE.set_differential( RHS )
  print( ODE.eqn_differential )
  ODE.set_initial( [ pymc.FFVar(1.) ] )
  ODE.set_quadrature( QUAD, Q )
  ODE.set_function( F )
  print( ODE.eqn_function )

  ODE.options.DISPLEVEL = 0
  ODE.options.INTMETH   = ODE.options.MSBDF
  ODE.options.NLINSOL   = ODE.options.NEWTON #FIXEDPOINT
  ODE.options.LINSOL    = ODE.options.DENSE  #DIAG
  ODE.setup()

  ODE.options.DISPLEVEL = 1
  stat = ODE.solve_state( [-1e0]*NS )

  return ODE

ODE = ode_define()
#ode_copy( ODE )
#ode_dag( ODE )
#ode_diff( ODE )
ode_dag2( ODE )

#ode2_define( 5 )

#help(cronos)
#help(cronos.ODESLV.solve_state)

