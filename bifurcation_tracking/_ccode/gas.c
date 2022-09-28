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
  const unsigned this_nodalind_c = 0;
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
    double this_intrp_d0t0_d1x0_c=0.0;
    double this_intrp_d0t0_d1x1_c=0.0;
    for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
    {
      this_intrp_d0t0_d1x0_c+= eleminfo->nodal_data[l_shape][this_nodalind_c][0] * shapeinfo->dx_shape_C2[ipt][l_shape][0];
      this_intrp_d0t0_d1x1_c+= eleminfo->nodal_data[l_shape][this_nodalind_c][0] * shapeinfo->dx_shape_C2[ipt][l_shape][1];
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
        BEGIN_RESIDUAL_CONTINUOUS_SPACE(eleminfo->nodal_local_eqn[l_test][this_nodalind_c],6.2831853071795862e+00*( dx_testfunction[l_test][0]*this_intrp_d0t0_d1x0_c+dx_testfunction[l_test][1]*this_intrp_d0t0_d1x1_c)*dx*this_intrp_d0t0_d0x_coordinate_x, shapeinfo->hanginfo_C2,this_nodalind_c,l_test)
          ADD_TO_RESIDUAL_CONTINUOUS_SPACE()
          BEGIN_JACOBIAN()
            for (unsigned int l_shape=0;l_shape<eleminfo->nnode_C2;l_shape++)
            {
              BEGIN_JACOBIAN_HANG(eleminfo->nodal_local_eqn[l_shape][this_nodalind_c], 6.2831853071795862e+00*( dx_testfunction[l_test][0]*shapeinfo->dx_shape_C2[ipt][l_shape][0]+dx_testfunction[l_test][1]*shapeinfo->dx_shape_C2[ipt][l_shape][1])*dx*this_intrp_d0t0_d0x_coordinate_x,shapeinfo->hanginfo_C2,this_nodalind_c,l_shape)
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
 functable->numfields_C2=1;
 functable->numfields_C2_bulk=1;
 functable->numfields_C2_basebulk=1;
 functable->fieldnames_C2=(char **)malloc(sizeof(char*)*functable->numfields_C2);
 SET_INTERNAL_FIELD_NAME(functable->fieldnames_C2,0, "c" );
 functable->dominant_space=strdup("C2");
 functable->max_dt_order=0;
 functable->moving_nodes=false;
 functable->num_res_jacs=1;
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

 functable->num_Z2_flux_terms = 0;
 functable->temporal_error_scales=calloc(7,sizeof(double)); 
 functable->discontinuous_refinement_exponents=calloc(7,sizeof(double));
 functable->num_ICs=1;
 functable->IC_names=(char**)calloc(functable->num_ICs,sizeof(char*));
 functable->InitialConditionFunc=(JITFuncSpec_InitialCondition_FiniteElement*)calloc(functable->num_ICs,sizeof(JITFuncSpec_InitialCondition_FiniteElement));
 SET_INTERNAL_FIELD_NAME(functable->IC_names,0, "" );
 functable->InitialConditionFunc[0]=&ElementalInitialConditions0;
 functable->DirichletConditionFunc=&ElementalDirichletConditions;
 functable->Dirichlet_set_size=4;
 functable->Dirichlet_set=(bool *)calloc(functable->Dirichlet_set_size,sizeof(bool)); 
 functable->Dirichlet_names=(char**)calloc(functable->Dirichlet_set_size,sizeof(char*));
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,0, "" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,1, "coordinate_y" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,2, "coordinate_x" );
 SET_INTERNAL_FIELD_NAME(functable->Dirichlet_names,3, "c" );
 functable->integration_order=0;
 functable->GeometricJacobian=&GeometricJacobian;
 my_func_table=functable;
}
