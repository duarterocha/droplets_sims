#define JIT_ELEMENT_SHARED_LIB
#include "jitbridge.h"

static JITFuncSpec_Table_FiniteElement_t * my_func_table;
#include "jitbridge_hang.h"


static void ResidualAndJacobian0(const JITElementInfo_t * eleminfo, const JITShapeInfo_t * shapeinfo,double * residuals, double *jacobian, double *mass_matrix,unsigned flag)
{
  int local_eqn, local_unknown;
  unsigned nummaster,nummaster2;
  double hang_weight,hang_weight2;
 double * t=shapeinfo->t;
 double * dt=shapeinfo->dt;
  const unsigned this_nodalind_velocity_x = 0;
  const unsigned this_nodalind_velocity_y = 1;
  const unsigned this_nodalind_pressure = 4;
  const unsigned this_nodalind_T = 2;
  const unsigned this_nodalind_streamfunc = 3;
  const unsigned this_nodalind_coordinate_x = 0;
  //START: Precalculate time derivatives of the necessary data
  PYOOMPH_AQUIRE_ARRAY(double, this_d1t0BDF2_T, eleminfo->nnode_C2)
  for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
  {
    this_d1t0BDF2_T[l_shape]=0.0;
    for (unsigned tindex=0;tindex<shapeinfo->timestepper_ntstorage;tindex++)
    {
      this_d1t0BDF2_T[l_shape] += shapeinfo->timestepper_weights_dt_BDF2_degr[tindex]*eleminfo->nodal_data[l_shape][this_nodalind_T][tindex];
    }
  }
  //END: Precalculate time derivatives of the necessary data

  //START: Spatial integration loop
  for(unsigned ipt=0;ipt<shapeinfo->n_int_pt;ipt++)
  {
    const double dx = shapeinfo->int_pt_weights[ipt];
    //START: Interpolate all required fields
    double this_intrp_d0t0_d0x_coordinate_x=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode;l_shape++)
    {
      this_intrp_d0t0_d0x_coordinate_x+= eleminfo->nodal_coords[l_shape][this_nodalind_coordinate_x][0] * shapeinfo->shape_Pos[ipt][l_shape];
    }
    double this_intrp_d0t0_d0x_velocity_x=0.0;
    double this_intrp_d0t0_d1x0_velocity_x=0.0;
    double this_intrp_d0t0_d1x1_velocity_x=0.0;
    double this_intrp_d0t0_d0x_velocity_y=0.0;
    double this_intrp_d0t0_d1x0_velocity_y=0.0;
    double this_intrp_d0t0_d1x1_velocity_y=0.0;
    double this_intrp_d0t0_d0x_T=0.0;
    double this_intrp_d0t0_d1x0_T=0.0;
    double this_intrp_d0t0_d1x1_T=0.0;
    double this_intrp_d1t0BDF2_d0x_T=0.0;
    double this_intrp_d0t0_d1x0_streamfunc=0.0;
    double this_intrp_d0t0_d1x1_streamfunc=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
    {
      this_intrp_d0t0_d0x_velocity_x+= eleminfo->nodal_data[l_shape][this_nodalind_velocity_x][0] * shapeinfo->shape_C2[ipt][l_shape];
      this_intrp_d0t0_d1x0_velocity_x+= eleminfo->nodal_data[l_shape][this_nodalind_velocity_x][0] * shapeinfo->dx_shape_C2[ipt][l_shape][0];
      this_intrp_d0t0_d1x1_velocity_x+= eleminfo->nodal_data[l_shape][this_nodalind_velocity_x][0] * shapeinfo->dx_shape_C2[ipt][l_shape][1];
      this_intrp_d0t0_d0x_velocity_y+= eleminfo->nodal_data[l_shape][this_nodalind_velocity_y][0] * shapeinfo->shape_C2[ipt][l_shape];
      this_intrp_d0t0_d1x0_velocity_y+= eleminfo->nodal_data[l_shape][this_nodalind_velocity_y][0] * shapeinfo->dx_shape_C2[ipt][l_shape][0];
      this_intrp_d0t0_d1x1_velocity_y+= eleminfo->nodal_data[l_shape][this_nodalind_velocity_y][0] * shapeinfo->dx_shape_C2[ipt][l_shape][1];
      this_intrp_d0t0_d0x_T+= eleminfo->nodal_data[l_shape][this_nodalind_T][0] * shapeinfo->shape_C2[ipt][l_shape];
      this_intrp_d0t0_d1x0_T+= eleminfo->nodal_data[l_shape][this_nodalind_T][0] * shapeinfo->dx_shape_C2[ipt][l_shape][0];
      this_intrp_d0t0_d1x1_T+= eleminfo->nodal_data[l_shape][this_nodalind_T][0] * shapeinfo->dx_shape_C2[ipt][l_shape][1];
      this_intrp_d1t0BDF2_d0x_T+= this_d1t0BDF2_T[l_shape] * shapeinfo->shape_C2[ipt][l_shape];
      this_intrp_d0t0_d1x0_streamfunc+= eleminfo->nodal_data[l_shape][this_nodalind_streamfunc][0] * shapeinfo->dx_shape_C2[ipt][l_shape][0];
      this_intrp_d0t0_d1x1_streamfunc+= eleminfo->nodal_data[l_shape][this_nodalind_streamfunc][0] * shapeinfo->dx_shape_C2[ipt][l_shape][1];
    }
    double this_intrp_d0t0_d0x_pressure=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C1;l_shape++)
    {
      this_intrp_d0t0_d0x_pressure+= eleminfo->nodal_data[l_shape][this_nodalind_pressure][0] * shapeinfo->shape_C1[ipt][l_shape];
    }
    //END: Interpolate all required fields


    // SUBEXPRESSIONS

 //Subexpressions // TODO: Check whether it is constant to take it out of the loop
    //Derivatives of subexpressions
    if (flag)
    {
  }
    //START: Contribution of the spaces
    double _res_contrib,_J_contrib;
    {
      double const * testfunction = shapeinfo->shape_C2[ipt];
      double * const * const dx_testfunction = shapeinfo->dx_shape_C2[ipt];
      double * const * const dX_testfunction = shapeinfo->dX_shape_C2[ipt];
      for (unsigned int l_test=0;l_test<eleminfo->nnode_C2;l_test++)
      {
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_T],2.0000000000000000e+00*dx*( 3.1415926535897931e+00*this_intrp_d0t0_d0x_coordinate_x*( this_intrp_d0t0_d1x0_T*dx_testfunction[l_test][0]+this_intrp_d0t0_d1x1_T*dx_testfunction[l_test][1])+3.1415926535897931e+00*( this_intrp_d0t0_d0x_velocity_x*this_intrp_d0t0_d1x0_T+this_intrp_d0t0_d0x_velocity_y*this_intrp_d0t0_d1x1_T)*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test]+3.1415926535897931e+00*this_intrp_d1t0BDF2_d0x_T*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test]), shapeinfo->hanginfo_C2,this_nodalind_T,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_x], 6.2831853071795862e+00*dx*shapeinfo->shape_C2[ipt][l_shape]*this_intrp_d0t0_d1x0_T*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test],shapeinfo->hanginfo_C2,this_nodalind_velocity_x,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_y], 6.2831853071795862e+00*dx*shapeinfo->shape_C2[ipt][l_shape]*this_intrp_d0t0_d1x1_T*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test],shapeinfo->hanginfo_C2,this_nodalind_velocity_y,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_T], 2.0000000000000000e+00*dx*( 3.1415926535897931e+00*this_intrp_d0t0_d0x_coordinate_x*( shapeinfo->dx_shape_C2[ipt][l_shape][0]*dx_testfunction[l_test][0]+shapeinfo->dx_shape_C2[ipt][l_shape][1]*dx_testfunction[l_test][1])+3.1415926535897931e+00*( this_intrp_d0t0_d0x_velocity_x*shapeinfo->dx_shape_C2[ipt][l_shape][0]+this_intrp_d0t0_d0x_velocity_y*shapeinfo->dx_shape_C2[ipt][l_shape][1])*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test]+3.1415926535897931e+00*shapeinfo->timestepper_weights_dt_BDF2_degr[0]*shapeinfo->shape_C2[ipt][l_shape]*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test]),shapeinfo->hanginfo_C2,this_nodalind_T,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
                ADD_TO_MASS_MATRIX_HANG_HANG(6.2831853071795862e+00*dx*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test])
              END_JACOBIAN_HANG()
            }
          END_JACOBIAN()
        END_RESIDUAL_CONTINUOUS_SPACE()
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_streamfunc],-2.0000000000000000e+00*dx*( 3.1415926535897931e+00*this_intrp_d0t0_d0x_coordinate_x*( this_intrp_d0t0_d1x0_streamfunc*dx_testfunction[l_test][0]+this_intrp_d0t0_d1x1_streamfunc*dx_testfunction[l_test][1])+3.1415926535897931e+00*( 2.0000000000000000e+00*this_intrp_d0t0_d0x_velocity_y-( this_intrp_d0t0_d1x1_velocity_x-this_intrp_d0t0_d1x0_velocity_y)*this_intrp_d0t0_d0x_coordinate_x)*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test]), shapeinfo->hanginfo_C2,this_nodalind_streamfunc,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_x], 6.2831853071795862e+00*dx*shapeinfo->dx_shape_C2[ipt][l_shape][1]*pow(this_intrp_d0t0_d0x_coordinate_x,2.0)*testfunction[l_test],shapeinfo->hanginfo_C2,this_nodalind_velocity_x,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_y], -6.2831853071795862e+00*( shapeinfo->dx_shape_C2[ipt][l_shape][0]*this_intrp_d0t0_d0x_coordinate_x+2.0000000000000000e+00*shapeinfo->shape_C2[ipt][l_shape])*dx*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test],shapeinfo->hanginfo_C2,this_nodalind_velocity_y,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_streamfunc], -6.2831853071795862e+00*dx*this_intrp_d0t0_d0x_coordinate_x*( shapeinfo->dx_shape_C2[ipt][l_shape][0]*dx_testfunction[l_test][0]+shapeinfo->dx_shape_C2[ipt][l_shape][1]*dx_testfunction[l_test][1]),shapeinfo->hanginfo_C2,this_nodalind_streamfunc,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
            }
          END_JACOBIAN()
        END_RESIDUAL_CONTINUOUS_SPACE()
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_velocity_x],6.2831853071795862e+00*dx*this_intrp_d0t0_d0x_coordinate_x*( dx_testfunction[l_test][0]*( 2.0000000000000000e+00*this_intrp_d0t0_d1x0_velocity_x-this_intrp_d0t0_d0x_pressure)+dx_testfunction[l_test][1]*( this_intrp_d0t0_d1x1_velocity_x+this_intrp_d0t0_d1x0_velocity_y)-( this_intrp_d0t0_d0x_pressure+-2.0000000000000000e+00*this_intrp_d0t0_d0x_velocity_x/this_intrp_d0t0_d0x_coordinate_x)/this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test]), shapeinfo->hanginfo_C2,this_nodalind_velocity_x,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_x], 6.2831853071795862e+00*dx*this_intrp_d0t0_d0x_coordinate_x*( 2.0000000000000000e+00*shapeinfo->shape_C2[ipt][l_shape]/pow(this_intrp_d0t0_d0x_coordinate_x,2.0)*testfunction[l_test]+2.0000000000000000e+00*shapeinfo->dx_shape_C2[ipt][l_shape][0]*dx_testfunction[l_test][0]+shapeinfo->dx_shape_C2[ipt][l_shape][1]*dx_testfunction[l_test][1]),shapeinfo->hanginfo_C2,this_nodalind_velocity_x,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_y], 6.2831853071795862e+00*dx*shapeinfo->dx_shape_C2[ipt][l_shape][0]*this_intrp_d0t0_d0x_coordinate_x*dx_testfunction[l_test][1],shapeinfo->hanginfo_C2,this_nodalind_velocity_y,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
            }
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C1;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_pressure], -6.2831853071795862e+00*dx*( shapeinfo->shape_C1[ipt][l_shape]/this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test]+shapeinfo->shape_C1[ipt][l_shape]*dx_testfunction[l_test][0])*this_intrp_d0t0_d0x_coordinate_x,shapeinfo->hanginfo_C1,this_nodalind_pressure,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
            }
          END_JACOBIAN()
        END_RESIDUAL_CONTINUOUS_SPACE()
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_velocity_y],-2.0000000000000000e+00*dx*( -3.1415926535897931e+00*this_intrp_d0t0_d0x_coordinate_x*( dx_testfunction[l_test][1]*( 2.0000000000000000e+00*this_intrp_d0t0_d1x1_velocity_y-this_intrp_d0t0_d0x_pressure)+dx_testfunction[l_test][0]*( this_intrp_d0t0_d1x1_velocity_x+this_intrp_d0t0_d1x0_velocity_y))+3.1415926535897931e+00*(*(my_func_table->global_parameters[0]))*this_intrp_d0t0_d0x_T*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test]), shapeinfo->hanginfo_C2,this_nodalind_velocity_y,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_x], 6.2831853071795862e+00*dx*shapeinfo->dx_shape_C2[ipt][l_shape][1]*this_intrp_d0t0_d0x_coordinate_x*dx_testfunction[l_test][0],shapeinfo->hanginfo_C2,this_nodalind_velocity_x,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_y], 6.2831853071795862e+00*dx*this_intrp_d0t0_d0x_coordinate_x*( shapeinfo->dx_shape_C2[ipt][l_shape][0]*dx_testfunction[l_test][0]+2.0000000000000000e+00*shapeinfo->dx_shape_C2[ipt][l_shape][1]*dx_testfunction[l_test][1]),shapeinfo->hanginfo_C2,this_nodalind_velocity_y,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_T], -6.2831853071795862e+00*dx*(*(my_func_table->global_parameters[0]))*shapeinfo->shape_C2[ipt][l_shape]*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test],shapeinfo->hanginfo_C2,this_nodalind_T,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
            }
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C1;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_pressure], -6.2831853071795862e+00*dx*shapeinfo->shape_C1[ipt][l_shape]*this_intrp_d0t0_d0x_coordinate_x*dx_testfunction[l_test][1],shapeinfo->hanginfo_C1,this_nodalind_pressure,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
            }
          END_JACOBIAN()
        END_RESIDUAL_CONTINUOUS_SPACE()
      }
    }
    {
      double const * testfunction = shapeinfo->shape_C1[ipt];
      double * const * const dx_testfunction = shapeinfo->dx_shape_C1[ipt];
      double * const * const dX_testfunction = shapeinfo->dX_shape_C1[ipt];
      for (unsigned int l_test=0;l_test<eleminfo->nnode_C1;l_test++)
      {
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_pressure],6.2831853071795862e+00*dx*this_intrp_d0t0_d0x_coordinate_x*( this_intrp_d0t0_d1x0_velocity_x+this_intrp_d0t0_d1x1_velocity_y+this_intrp_d0t0_d0x_velocity_x/this_intrp_d0t0_d0x_coordinate_x)*testfunction[l_test], shapeinfo->hanginfo_C1,this_nodalind_pressure,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_x], 6.2831853071795862e+00*( shapeinfo->dx_shape_C2[ipt][l_shape][0]+shapeinfo->shape_C2[ipt][l_shape]/this_intrp_d0t0_d0x_coordinate_x)*dx*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test],shapeinfo->hanginfo_C2,this_nodalind_velocity_x,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_velocity_y], 6.2831853071795862e+00*dx*shapeinfo->dx_shape_C2[ipt][l_shape][1]*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test],shapeinfo->hanginfo_C2,this_nodalind_velocity_y,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
            }
          END_JACOBIAN()
        END_RESIDUAL_CONTINUOUS_SPACE()
      }
    }
    //END: Contribution of the spaces
  }
  //END: Spatial integration loop

}



//Derivative wrt. global parameter <global param: Ra>
static void dResidual0dParameter_0(const JITElementInfo_t * eleminfo, const JITShapeInfo_t * shapeinfo,double * residuals, double *jacobian, double *mass_matrix,unsigned flag)
{
  int local_eqn, local_unknown;
  unsigned nummaster,nummaster2;
  double hang_weight,hang_weight2;
 double * t=shapeinfo->t;
 double * dt=shapeinfo->dt;
  const unsigned this_nodalind_velocity_y = 1;
  const unsigned this_nodalind_T = 2;
  const unsigned this_nodalind_coordinate_x = 0;
  //START: Precalculate time derivatives of the necessary data
  //END: Precalculate time derivatives of the necessary data

  //START: Spatial integration loop
  for(unsigned ipt=0;ipt<shapeinfo->n_int_pt;ipt++)
  {
    const double dx = shapeinfo->int_pt_weights[ipt];
    //START: Interpolate all required fields
    double this_intrp_d0t0_d0x_coordinate_x=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode;l_shape++)
    {
      this_intrp_d0t0_d0x_coordinate_x+= eleminfo->nodal_coords[l_shape][this_nodalind_coordinate_x][0] * shapeinfo->shape_Pos[ipt][l_shape];
    }
    double this_intrp_d0t0_d0x_T=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
    {
      this_intrp_d0t0_d0x_T+= eleminfo->nodal_data[l_shape][this_nodalind_T][0] * shapeinfo->shape_C2[ipt][l_shape];
    }
    //END: Interpolate all required fields


    // SUBEXPRESSIONS

 //Subexpressions // TODO: Check whether it is constant to take it out of the loop
    //Derivatives of subexpressions
    if (flag)
    {
  }
    //START: Contribution of the spaces
    double _res_contrib,_J_contrib;
    {
      double const * testfunction = shapeinfo->shape_C2[ipt];
      double * const * const dx_testfunction = shapeinfo->dx_shape_C2[ipt];
      double * const * const dX_testfunction = shapeinfo->dX_shape_C2[ipt];
      for (unsigned int l_test=0;l_test<eleminfo->nnode_C2;l_test++)
      {
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_velocity_y],-6.2831853071795862e+00*dx*this_intrp_d0t0_d0x_T*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test], shapeinfo->hanginfo_C2,this_nodalind_velocity_y,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_T], -6.2831853071795862e+00*dx*shapeinfo->shape_C2[ipt][l_shape]*this_intrp_d0t0_d0x_coordinate_x*testfunction[l_test],shapeinfo->hanginfo_C2,this_nodalind_T,l_shape)
                ADD_TO_JACOBIAN_HANG_HANG()
              END_JACOBIAN_HANG()
            }
          END_JACOBIAN()
        END_RESIDUAL_CONTINUOUS_SPACE()
      }
    }
    //END: Contribution of the spaces
  }
  //END: Spatial integration loop

}

// INITIAL CONDITION 
static double ElementalInitialConditions0(const JITElementInfo_t * eleminfo, int field_index,double *_x, double *_xlagr,double t,int flag,double default_val)
{
  return default_val;
}

static double ElementalDirichletConditions(const JITElementInfo_t * eleminfo, int field_index,double *_x, double *_xlagr,double t,double default_val)
{
  return default_val;
}

static double GeometricJacobian(const JITElementInfo_t * eleminfo, const double * _x)
{
  return 6.2831853071795862e+00*_x[0];
}

static void GetZ2Fluxes(const JITElementInfo_t * eleminfo, const JITShapeInfo_t * shapeinfo, double * Z2Flux)
{
 const unsigned ipt=0;  //The shape information is evaled at the zero'th integration point storage, which is filled with the local info before

  const unsigned this_nodalind_streamfunc = 3;
  //START: Precalculate time derivatives of the necessary data
  //END: Precalculate time derivatives of the necessary data

    //START: Interpolate all required fields
    double this_intrp_d0t0_d1x0_streamfunc=0.0;
    double this_intrp_d0t0_d1x1_streamfunc=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
    {
      this_intrp_d0t0_d1x0_streamfunc+= eleminfo->nodal_data[l_shape][this_nodalind_streamfunc][0] * shapeinfo->dx_shape_C2[ipt][l_shape][0];
      this_intrp_d0t0_d1x1_streamfunc+= eleminfo->nodal_data[l_shape][this_nodalind_streamfunc][0] * shapeinfo->dx_shape_C2[ipt][l_shape][1];
    }
    //END: Interpolate all required fields

  Z2Flux[0] = this_intrp_d0t0_d1x0_streamfunc;
  Z2Flux[1] = this_intrp_d0t0_d1x1_streamfunc;
}


static double EvalIntegralExpression(const JITElementInfo_t * eleminfo, const JITShapeInfo_t * shapeinfo, unsigned index)
{
  double * t=shapeinfo->t;
  double * dt=shapeinfo->dt;

  const unsigned this_nodalind_velocity_x = 0;
  const unsigned this_nodalind_velocity_y = 1;
  const unsigned this_nodalind_streamfunc = 3;
  const unsigned this_nodalind_coordinate_x = 0;
  //START: Precalculate time derivatives of the necessary data
  //END: Precalculate time derivatives of the necessary data

  double res=0.0;
  for(unsigned ipt=0;ipt<shapeinfo->n_int_pt;ipt++)
  {
    //START: Interpolate all required fields
    double this_intrp_d0t0_d0x_coordinate_x=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode;l_shape++)
    {
      this_intrp_d0t0_d0x_coordinate_x+= eleminfo->nodal_coords[l_shape][this_nodalind_coordinate_x][0] * shapeinfo->shape_Pos[ipt][l_shape];
    }
    double this_intrp_d0t0_d0x_velocity_x=0.0;
    double this_intrp_d0t0_d0x_velocity_y=0.0;
    double this_intrp_d0t0_d0x_streamfunc=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
    {
      this_intrp_d0t0_d0x_velocity_x+= eleminfo->nodal_data[l_shape][this_nodalind_velocity_x][0] * shapeinfo->shape_C2[ipt][l_shape];
      this_intrp_d0t0_d0x_velocity_y+= eleminfo->nodal_data[l_shape][this_nodalind_velocity_y][0] * shapeinfo->shape_C2[ipt][l_shape];
      this_intrp_d0t0_d0x_streamfunc+= eleminfo->nodal_data[l_shape][this_nodalind_streamfunc][0] * shapeinfo->shape_C2[ipt][l_shape];
    }
    //END: Interpolate all required fields

    const double dx = shapeinfo->int_pt_weights[ipt];
    switch (index)
    {
      case 0 :  res+= dx; break; // _cross_area [ 1 ]
      case 1 :  res+= dx*step(this_intrp_d0t0_d0x_streamfunc); break; // _pos_stream_area [ 1 ]
      case 2 :  res+= dx*this_intrp_d0t0_d0x_coordinate_x*( pow(this_intrp_d0t0_d0x_velocity_x,2.0)+pow(this_intrp_d0t0_d0x_velocity_y,2.0)); break; // _v_sqr_int [ 2*Pi ]
      case 3 :  res+= dx*this_intrp_d0t0_d0x_coordinate_x; break; // _volume [ 2*Pi ]
    }
  }
   return res;
}


JIT_API void JIT_ELEMENT_init(JITFuncSpec_Table_FiniteElement_t *functable)
{
 functable->nodal_dim=2;
 functable->lagr_dim=2;
 functable->fd_jacobian=false; 
 functable->fd_position_jacobian=false; 
 functable->debug_jacobian_epsilon = 0;
 functable->stop_on_jacobian_difference = false;
 functable->numfields_Pos=4;
 functable->fieldnames_Pos=(char **)malloc(sizeof(char*)*functable->numfields_Pos);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,0, "coordinate_x" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,1, "coordinate_y" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,2, "lagrangian_x" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_Pos,3, "lagrangian_y" );
 functable->numfields_C2=4;
 functable->numfields_C2_bulk=4;
 functable->numfields_C2_basebulk=4;
 functable->fieldnames_C2=(char **)malloc(sizeof(char*)*functable->numfields_C2);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,0, "velocity_x" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,1, "velocity_y" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,2, "T" );
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,3, "streamfunc" );
 functable->numfields_C1=1;
 functable->numfields_C1_bulk=1;
 functable->numfields_C1_basebulk=1;
 functable->fieldnames_C1=(char **)malloc(sizeof(char*)*functable->numfields_C1);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C1,0, "pressure" );
 functable->dominant_space=strdup("C2");
 functable->max_dt_order=1;
 functable->moving_nodes=false;
 functable->num_res_jacs=1;
 functable->numglobal_params=1;
 functable->global_paramindices=(unsigned *)malloc(sizeof(unsigned)*functable->numglobal_params);
 functable->global_parameters=(double **)calloc(functable->numglobal_params,sizeof(double*));
 functable->global_paramindices[0]=0;
 functable->ParameterDerivative=(JITFuncSpec_ResidualAndJacobian_FiniteElement **)malloc(sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement)*functable->num_res_jacs);
 functable->ParameterDerivative[0]=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)malloc(sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement)*functable->numglobal_params);
 functable->ParameterDerivative[0][0]=&dResidual0dParameter_0;
 functable->ResidualAndJacobian_NoHang=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->ResidualAndJacobian=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->ResidualAndJacobianSteady=(JITFuncSpec_ResidualAndJacobian_FiniteElement *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_ResidualAndJacobian_FiniteElement));
 functable->shapes_required_ResJac=(JITFuncSpec_RequiredShapes_FiniteElement_t *)calloc(functable->num_res_jacs,sizeof(JITFuncSpec_RequiredShapes_FiniteElement_t));
 functable->res_jac_names=(char**)calloc(functable->num_res_jacs,sizeof(char*));
 SET_INTERNAL_FIELD_NAME(functable->res_jac_names,0, "" );
 functable->ResidualAndJacobian_NoHang[0]=&ResidualAndJacobian0;
 functable->ResidualAndJacobian[0]=&ResidualAndJacobian0;
 functable->ResidualAndJacobianSteady[0]=&ResidualAndJacobian0;
  functable->shapes_required_ResJac[0].psi_Pos = true;
  functable->shapes_required_ResJac[0].dx_psi_C2 = true;
  functable->shapes_required_ResJac[0].psi_C2 = true;
  functable->shapes_required_ResJac[0].psi_C1 = true;

 functable->num_Z2_flux_terms = 2;
 functable->GetZ2Fluxes=&GetZ2Fluxes;
  functable->shapes_required_Z2Fluxes.dx_psi_C2 = true;
 functable->temporal_error_scales=calloc(11,sizeof(double)); 
 functable->discontinuous_refinement_exponents=calloc(11,sizeof(double));
 functable->num_ICs=1;
 functable->IC_names=(char**)calloc(functable->num_ICs,sizeof(char*));
 functable->InitialConditionFunc=(JITFuncSpec_InitialCondition_FiniteElement*)calloc(functable->num_ICs,sizeof(JITFuncSpec_InitialCondition_FiniteElement));
 SET_INTERNAL_FIELD_NAME(functable->IC_names,0, "" );
 functable->InitialConditionFunc[0]=&ElementalInitialConditions0;
 functable->DirichletConditionFunc=&ElementalDirichletConditions;
 functable->Dirichlet_set_size=8;
 functable->Dirichlet_set=(bool *)calloc(functable->Dirichlet_set_size,sizeof(bool)); 
 functable->Dirichlet_names=(char**)calloc(functable->Dirichlet_set_size,sizeof(char*));
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,0, "" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,1, "coordinate_y" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,2, "coordinate_x" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,3, "velocity_x" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,4, "velocity_y" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,5, "T" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,6, "streamfunc" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,7, "pressure" );
 functable->numintegral_expressions=4;
 functable->integral_expressions_names=(char **)malloc(sizeof(char*)*functable->numintegral_expressions);
 SET_INTERNAL_FIELD_NAME(functable->integral_expressions_names,0,"_cross_area");
 SET_INTERNAL_FIELD_NAME(functable->integral_expressions_names,1,"_pos_stream_area");
 SET_INTERNAL_FIELD_NAME(functable->integral_expressions_names,2,"_v_sqr_int");
 SET_INTERNAL_FIELD_NAME(functable->integral_expressions_names,3,"_volume");
 functable->EvalIntegralExpression=&EvalIntegralExpression;
  functable->shapes_required_IntegralExprs.psi_Pos = true;
  functable->shapes_required_IntegralExprs.psi_C2 = true;
 functable->integration_order=0;
 functable->GeometricJacobian=&GeometricJacobian;
 my_func_table=functable;
}
