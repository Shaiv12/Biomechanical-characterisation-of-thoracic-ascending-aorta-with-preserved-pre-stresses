%% ckOptimize_StraightTube_2stepOptimize


clear; close all; clc;
tic
fontSize=30;
faceAlpha1=0.8;
% markerSize=40;
markerSize2=20;
lineWidth=3;

% Path names

savePath='C:\Users\p70068216\Desktop\ckOptimize_straightTube';

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt'];%Log file name for exporting stress

%Specifying geometry parameters
pointSpacing= 3;
lenParam = 6;
radiusOuter1=15;
% radiusInner1=9;
radiusOuter2=15;
% radiusInner2=7;
vesselLength= 2*lenParam*radiusOuter1;
% vesselLength= 100;
wallThickness = 1;
numThickenSteps = 1;



%Load
hG2N = 0.000133322;
pressureDiasmmHg = 80;
pressureSysmmHg = 130;
appliedPressureDias= (pressureDiasmmHg*hG2N);
appliedPressureSys= (pressureSysmmHg*hG2N);



%Material parameter set
c1=0.9; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=12; %Optimum number of iterations
max_retries=10; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

nRad=ceil((2*pi*mean([radiusOuter1 radiusOuter2]))/pointSpacing); %Number of radial steps

nLong = ceil(vesselLength/pointSpacing);

nLong = nLong + iseven(nLong);

% Creating input structure
inputStruct.cylRadius= radiusOuter1;
inputStruct.numRadial=nRad;
inputStruct.cylHeight=vesselLength;
inputStruct.numHeight=nLong;
inputStruct.meshType='quad';

% Derive patch data for a cylinder
[F,V]=patchcylinder(inputStruct);


C = ones(size(F,1),1);
cFigure;hold on;
gpatch(F,V,C,'k',0.85)
axisGeom;
gdrawnow;


%% Inverting offset direction
%Converting face thickness data to vertex data

[~,~,N]=patchNormal(F,V);

%%
figStruct.ColorDef='white';
figStruct.Color='w';

cFigure(figStruct); hold on;
gpatch(F,V,C,'kw',0.5);
title('Mesh Offset')
quiverVec(V,-N.*wallThickness,[],wallThickness);
axisGeom;


%% Thicken to create hexahedral elements
[ET,VT,Fp1,Fp2]=patchThick(F,V,-1,wallThickness,numThickenSteps);

CT=repmat(C,numThickenSteps,1);

ET=ET(:,[5:8 1:4]);
FT_inner=Fp2;
FT_outer=Fp1;
FT_inner = fliplr(FT_inner);
indInner = unique(FT_inner);
indOuter =unique(FT_outer);
indRingAsc = indInner(VT(indInner,3)<eps(0) & VT(indInner,3)> -eps(0));
indRingAsc_outer = indOuter(VT(indOuter,3)<eps(0) & VT(indOuter,3)> -eps(0));

%%
cFigure(figStruct); hold on;
gpatch(FT_inner,VT,'w','kw',0.5);
plotV(VT(indRingAsc,:),'r.','Markersize',25);
axisGeom;
gdrawnow;
%%


F_pressure = FT_inner;

[~,logicPositive]=hexVol(ET,VT);
% logicPositive
if any(logicPositive==0)
    error('Negative hex volume found');
end


[FT,CFT,CFT_type]=element2patch(ET,CT);
indBoundary = tesBoundary(FT,VT);
FT_b= FT(indBoundary,:);
CFT_type_b=CFT_type(indBoundary,:);
CFT_b= CFT(indBoundary,:);
bcSupportList=unique(FT_b(ismember(CFT_type_b,[3 4]),:));%list of indices for all boundary nodes (whole geometry)

bcnodMemb0= ismember(VT,VT(bcSupportList,:),'rows'); %filter out nodes common to only boundaries and feature 1
VT0 = VT(bcnodMemb0,:);
cFigure(figStruct); hold on;
gpatch(FT_b,VT,CFT_type_b,'k');
plotV(VT0,'k.','Markersize',25);
% gpatch(FT,VT,CFT,'k');
axisGeom;
camlight headlight;
icolorbar;
gdrawnow;



%% Acquiring diastolic and systolic nodal displacmenet wrt undeformed state

appliedPressure = [appliedPressureDias appliedPressureSys];

for q = 1:1:size(appliedPressure,2)
    
    %% FEBio input struct
    
    %Get a template with default settings
    [febio_spec]=febioStructTemplate;
    
    %febio_spec version
    febio_spec.ATTR.version='3.0';
    
    %Module section
    febio_spec.Module.ATTR.type='solid';
    
    %Control section
    febio_spec.Control.analysis='STATIC';
    febio_spec.Control.time_steps=numTimeSteps;
    febio_spec.Control.step_size=1/numTimeSteps;
    febio_spec.Control.solver.max_refs=max_refs;
    febio_spec.Control.solver.max_ups=max_ups;
    febio_spec.Control.time_stepper.dtmin=dtmin;
    febio_spec.Control.time_stepper.dtmax=dtmax;
    febio_spec.Control.time_stepper.max_retries=max_retries;
    febio_spec.Control.time_stepper.opt_iter=opt_iter;
    
    %Material section
    materialName1='Material1';
    febio_spec.Material.material{1}.ATTR.name=materialName1;
    febio_spec.Material.material{1}.ATTR.type='Ogden';
    febio_spec.Material.material{1}.ATTR.id=1;
    febio_spec.Material.material{1}.c1=c1;
    febio_spec.Material.material{1}.m1=m1;
    febio_spec.Material.material{1}.c2=c1;
    febio_spec.Material.material{1}.m2=-m1;
    febio_spec.Material.material{1}.k=k;
    
    %Mesh section
    % -> Nodes
    febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(VT,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL=VT; %The nodel coordinates
    
    % -> Elements
    partName1='Part1';
    febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
    febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type
    febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(ET,1))'; %Element id's
    febio_spec.Mesh.Elements{1}.elem.VAL=ET; %The element matrix
    
    % -> Surfaces
    surfaceName1='LoadedSurface';
    febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
    febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
    febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure;
    
    % -> NodeSets
    nodeSetName1='bcSupportList';
    febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
    febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);
    
    
    %MeshDomains section
    febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
    febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;
    
    %Boundary condition section
    % -> Fix boundary conditions
    febio_spec.Boundary.bc{1}.ATTR.type='fix';
    febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{1}.dofs='x,y,z';
    
    %Loads section
    % -> Surface load
    % -> Surface load
    febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
    febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
    febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
    % febio_spec.Loads.surface_load{1}.pressure.VAL=appliedPressureSys;
    febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;
    febio_spec.Loads.surface_load{1}.pressure.VAL=1;
    
    %LoadData section
    % -> load_controller
    febio_spec.LoadData.load_controller{1}.ATTR.id=1;
    febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
    febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
    febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 appliedPressure(q)];
    
    
    %Output section
    % -> log file
    febio_spec.Output.logfile.ATTR.file=febioLogFileName;
    febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
    febio_spec.Output.logfile.node_data{1}.VAL=1:size(VT,1);
    
    
    febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
    febio_spec.Output.logfile.element_data{1}.ATTR.data='sx;sy;sz';
    febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
    febio_spec.Output.logfile.element_data{1}.VAL=1:size(ET,1);
    
    
    %% Exporting the FEBio input file
    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
    
    %% Running FEBio
    febioAnalysis.run_filename=febioFebFileName; %The input file name
    febioAnalysis.run_logname=febioLogFileName; %The name for the log file
    febioAnalysis.disp_on=1; %Display information on the command window
    febioAnalysis.disp_log_on=1; %Display convergence information in the command window
    febioAnalysis.runMode='external';%'internal';
    febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
    febioAnalysis.maxtpi=1e99; %Max analysis time
    febioAnalysis.maxLogCheckTime=10; %Max log file checking time
    
    [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!
    
    
    
    %% Import FEBio results
    
    if runFlag==1 %i.e. a succesful run
        
        %Importing displacement
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
        N_disp_mat=dataStruct.data; %Displacement
        U{q} = N_disp_mat(:,:,end);
        
        
        
        %import elemental stresses from logfile
        dataStruct_stress = importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
        E_stress_mat{q}= dataStruct_stress.data; %Elemental stresses
    end
end






%% Creating diastolic and systolic geomnetries

V_DEF_dias = VT + U{1,1};
V_DEF_sys = VT + U{1,2};

meanDiasSynthStress = mean(sqrt(sum(E_stress_mat{1, 1}(:,:,end).^2,2)));


% %%
% cFigure;hold on;
% subplot(1,3,1);hold on;
% gpatch(FT,VT,'rw','k',1);
% % gpatch(FT_inner,VT,'rw','none',1);
% % plotV(VT(indRingAsc,:),'g-','Linewidth',3);
% plotV(VT,'g.','Markersize',20)
% plotV(VT0,'k.','Markersize',20);
% axis off
% axisGeom
% 
% subplot(1,3,2);hold on;
% gpatch(FT,VT,'rw','none',1);
% % gpatch(FT,VT,'rw','k',1);
% gpatch(FT,V_DEF_dias,'none','k',1);
% % plotV(V_DEF_dias(indRingAsc,:),'g-','Linewidth',3);
% plotV(V_DEF_dias,'g.','Markersize',20)
% plotV(VT0,'k.','Markersize',20);
% axis off
% axisGeom
% 
% 
% 
% subplot(1,3,3);hold on;
% gpatch(FT,VT,'rw','none',1);
% % gpatch(FT,VT,'rw','k',1);
% gpatch(FT,V_DEF_sys,'none','k',1);
% % plotV(V_DEF_dias(indRingAsc,:),'g-','Linewidth',3);
% plotV(V_DEF_sys,'g.','Markersize',20)
% plotV(VT0,'k.','Markersize',20);
% axis off
% axisGeom

%%

% subplot(1,3,3);hold on;
% gpatch(FT,V_DEF_sys,'none','k',1);
% gpatch(FT_inner,V_DEF_sys,'rw','none',1);
% plotV(V_DEF_sys(indRingAsc,:),'g-','Linewidth',3);
% plotV(VT0,'k.','Markersize',20);
% axis off
% axisGeom

% cFigure; hold on;
% gpatch(FT,V_DEF_sys,'none','k',1);
% gpatch(FT,V_DEF_dias,'rw','none',0.3);
% % plotV(V_DEF_dias,'g.','Markersize',15);
% plotV(V_DEF_sys,'g.','Markersize',15);
% V_def = V_DEF_sys-V_DEF_dias;
% quiverVec(V_DEF_dias,V_def,norm((V_def(2,:))/(norm(V_def(2,:)))),'b')
% axis off
% axisGeom

%%


%% Create elements required for optimisation
VT_img = V_DEF_dias; %Diastolic geomtery

%Diastolic middle ring
vDiasInnerRingNodes = V_DEF_dias(indRingAsc,:);%Inner ring
vDiasOuterRingNodes = V_DEF_dias(indRingAsc_outer,:);%outer ring

meanDiasAsc = mean(vDiasInnerRingNodes);
vecDiasAsc = vDiasInnerRingNodes - meanDiasAsc;
r_exp_Dias = mean(sqrt(sum(vecDiasAsc.^2,2))); %experimental averaged inner radius of the tube in diastole

%Systolic middle ring
vSysInnerRingNodes = V_DEF_sys(indRingAsc,:);
vSysOuterRingNodes = V_DEF_sys(indRingAsc_outer,:);
meanSysAsc = mean(vSysInnerRingNodes);
vecSysAsc = vSysInnerRingNodes - meanSysAsc;
r_exp_Sys = mean(sqrt(sum(vecSysAsc.^2,2))); %experimental averaged inner radius of the tube in systole


%Operations on outer radius (for thickness)
meanDiasAsc_outer = mean(vDiasOuterRingNodes);
vecDiasAsc_outer = vDiasOuterRingNodes - meanDiasAsc_outer;
r_exp_Dias_outer = mean(sqrt(sum(vecDiasAsc_outer.^2,2))); %experimental averaged outer radius of the tube in diastole

meanSysAsc_outer = mean(vSysOuterRingNodes);
vecSysAsc_outer = vSysOuterRingNodes - meanSysAsc_outer;
r_exp_Sys_outer = mean(sqrt(sum(vecSysAsc_outer.^2,2))); %experimental averaged outer radius



%Thickness values
hDias = r_exp_Dias_outer-r_exp_Dias;
hSys = r_exp_Sys_outer- r_exp_Sys;


%Stress value in systolic configuration (Laplace's law)

% sigmatt_Dias = (appliedPressureDias*r_exp_Dias)/hDias;
% sigmatt_Sys = (appliedPressureSys*r_exp_Sys)/hSys;

stress_expDias = max(max(E_stress_mat{1,1}(:,:,end)));
stress_expSys = max(max(E_stress_mat{1,2}(:,:,end)));





%% Acquire nodal displacements between diastolic and systolic geometry
U_find = V_DEF_sys-V_DEF_dias;


%% Creating Objetcive struct for getting optimization and determining U_find

objectiveStruct.r_exp_Sys = r_exp_Sys;
objectiveStruct.r_exp_Dias = r_exp_Dias;
objectiveStruct.innerRingInd = indRingAsc;
objectiveStruct.VT_img = VT_img;

objectiveStruct.appliedPressureSys = appliedPressureSys;
objectiveStruct.appliedPressureDias = appliedPressureDias;

objectiveStruct.febio_spec = febio_spec;
objectiveStruct.febioAnalysis= febioAnalysis;
objectiveStruct.febioFebFileName = febioFebFileName;
objectiveStruct.savePath = savePath;
objectiveStruct.febioLogFileName_disp= febioLogFileName_disp;
objectiveStruct.febioLogFileName_stress=febioLogFileName_stress;

objectiveStruct.U = U_find; %Add displacement field element to objective structure for optimization
objectiveStruct.laplaceStress= stress_expSys;
objectiveStruct.laplaceStress_dias=stress_expDias;




%%

%Optimisation settings
maxNumberIterations=1000; %Maximum number of optimization iterations
maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=1e-6; %Tolerance on objective function value
parameterTolerance=1e-6; %Tolerance on parameter variation
displayTypeIterations='iter';


SPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
SPT_options = optimoptions(SPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
    'MaxIter',maxNumberIterations,...
    'Display',displayTypeIterations,...
    'FinDiffRelStep',1e-2,...
    'DiffMaxChange',0.5);



OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
    'MaxIter',maxNumberIterations,...
    'Display',displayTypeIterations,...
    'FinDiffRelStep',1e-2,...
    'FinDiffType','central',...
    'DiffMaxChange',0.5);



%% Start parameter loop here
counter = 1;
Cinit= 0.5;
kinit = 1;
res= 1;

while res>0.01
    
    
    switch counter
        case 1
            S= kinit;
            C= Cinit;
            
        case counter>1
          S=K_new;
          C=C_new;
    end
    
    
    
    %% Creating undeformed state 
    
    
    Pn = C;
    [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveBackwardInc(Pn,objectiveStruct,r_exp_Sys,r_exp_Dias,S),Pn,[],[],OPT_options);
    
    C= Pn_opt;
    
    %Get shrinkage factor for a given material parameter
      Sn= S;
%     Sn  = 1;
    [Sn_opt,SPT_out.resnorm,SPT_out.residual]= lsqnonlin(@(Sn) getShrinkFact(Sn,objectiveStruct,C),Sn,[],[],SPT_options);
    
    S= Sn_opt;
    
    %% Creating undeformed state
    V_unDEF = VT_img - S*U_find;
    
    indVer=indRingAsc;
    V_unDEF_innerRing = V_unDEF(indVer,:);
    mean_V_unDEF = mean(V_unDEF_innerRing);
    allRadialVecs_unDEF = V_unDEF_innerRing - (mean_V_unDEF).*ones (size(V_unDEF_innerRing));
    meanRadius_unDEF = mean(sqrt(sum( allRadialVecs_unDEF.^2,2)));
    
    %% Run simulations in two steps
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_unDEF,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL=V_unDEF; %The nodel coordinates
    febio_spec.Output.logfile.node_data{1}.VAL=1:size(V_unDEF,1);
    
    febio_spec.Material.material{1}.c1= C; %setting new stiffness
    
    %Inflate till diastole
    febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureDias];
    
    %Exporting the FEBio input file
    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
    
    %START FEBio
    [runFlag]=runMonitorFEBio(febioAnalysis);
    
    
    if runFlag==1
        
        %Importing displacement
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
        N_disp_mat=dataStruct.data; %Displacement
        U_dias = N_disp_mat(:,:,end);
        %import elemental stresses from logfile
        dataStruct_stress = importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
        E_stress_mat= dataStruct_stress.data; %Elemental stresses
        
        %     stress_simDias = max(max(E_stress_mat(:,:,end)));
        
        %Derive deformed configuration
        V_DEF_dias= V_unDEF + U_dias;
        
        
        %Obtained systolic outer ring and mean radius
        V_DEF_dias_innerRing = V_DEF_dias(indVer,:);
        mean_V_DEF_dias = mean(V_DEF_dias_innerRing);
        allRadialVecs_dias = V_DEF_dias_innerRing - (mean_V_DEF_dias).*ones (size(V_DEF_dias_innerRing));
        meanRadius_dias = mean(sqrt(sum( allRadialVecs_dias.^2,2)));
        
    end
    
    
    %Inflate till systole
    febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureSys];
    
    %Exporting the FEBio input file
    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
    
    %START FEBio
    [runFlag]=runMonitorFEBio(febioAnalysis);
    
    
    if runFlag==1
        
        %Importing displacement
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
        N_disp_mat=dataStruct.data; %Displacement
        U_sys = N_disp_mat(:,:,end);
        V_DEF_sys= V_unDEF + U_sys;
        
        
        %Obtained systolic outer ring and mean radius
        V_DEF_sys_innerRing = V_DEF_sys(indVer,:);
        mean_V_DEF_sys = mean(V_DEF_sys_innerRing);
        allRadialVecs_sys = V_DEF_sys_innerRing - (mean_V_DEF_sys).*ones (size(V_DEF_sys_innerRing));
        meanRadius_sys = mean(sqrt(sum( allRadialVecs_sys.^2,2)));
    end
    
    
    K_new=S;
    C_new= C;
    
    %     res = abs((r_exp_Sys - meanRadius_sys ) + (r_exp_Dias-meanRadius_dias))
    res = abs((r_exp_Sys - r_exp_Dias ) - (meanRadius_sys-meanRadius_dias))
    counter= counter+1;
end


save straightTubeOptimize15.mat  


%% FUNCTIONS


function [Sopt,SPT_stats_out]=getShrinkFact(Sn,objectiveStruct,C)

%% Get variables
r_exp_Dias=objectiveStruct.r_exp_Dias;
indVer=objectiveStruct.innerRingInd;
VT_img = objectiveStruct.VT_img;

appliedPressureDias = objectiveStruct.appliedPressureDias ;

febio_spec=objectiveStruct.febio_spec ;
febioAnalysis= objectiveStruct.febioAnalysis;
febioFebFileName= objectiveStruct.febioFebFileName;
savePath= objectiveStruct.savePath;
febioLogFileName_disp = objectiveStruct.febioLogFileName_disp;
febioLogFileName_stress=objectiveStruct.febioLogFileName_stress;

% P= Pn.*objectiveStruct.parNormFactors;
U_find= objectiveStruct.U;
stress_expDias =objectiveStruct.laplaceStress_dias;

%% Get optimization parameters
K = Sn;

%% Creating undeformed state
V_unDEF = VT_img - K*U_find;

%%
% disp('SETTING PARAMETERS...');
% disp(['Proposed (norm.): ',sprintf(repmat('%6.16e ',[1,numel(Pn)]),Pn)]);
% disp(['Set (constr.)   : ',sprintf(repmat('%6.16e ',[1,numel(P)]),P)]);

%% Inflate till systole
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_unDEF,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V_unDEF; %The nodel coordinates
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V_unDEF,1);

febio_spec.Material.material{1}.c1= C; %setting new stiffness
febio_spec.Material.material{1}.c2= C;
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureDias];

%Exporting the FEBio input file
febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%START FEBio
[runFlag]=runMonitorFEBio(febioAnalysis);


if runFlag==1
    
    %Importing displacement
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    N_disp_mat=dataStruct.data; %Displacement
    U_dias = N_disp_mat(:,:,end);
    %import elemental stresses from logfile
    dataStruct_stress = importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
    E_stress_mat= dataStruct_stress.data; %Elemental stresses
    
    stress_simDias = max(max(E_stress_mat(:,:,end)));
    
    %Derive deformed configuration
    V_DEF_dias= V_unDEF + U_dias;
    
    
    
    
    %Obtained systolic outer ring and mean radius
    V_DEF_dias_innerRing = V_DEF_dias(indVer,:);
    mean_V_DEF_sys = mean(V_DEF_dias_innerRing);
    allRadialVecs_sys = V_DEF_dias_innerRing - (mean_V_DEF_sys).*ones (size(V_DEF_dias_innerRing));
    meanRadius_dias = mean(sqrt(sum( allRadialVecs_sys.^2,2)));
    
    
    %         Sopt = abs(r_exp_Dias - meanRadius_dias) + abs(stress_simDias-stress_expDias);
    
    
    vDev = VT_img - V_DEF_dias;
    Sopt = (sqrt(sum(vDev.^2,2)));
    %     Sopt = vDev;
    
    
    
    SPT_stats_out.Sopt=Sopt;
    
    
    
    
else
    Sopt=NaN(2,1); %Squared differences
    SPT_stats_out=[];
end
end




function [Popt,OPT_stats_out]=objectiveBackwardInc(Pn,objectiveStruct,r_exp_Sys,r_exp_Dias,S)

%% Get variables
indVer=objectiveStruct.innerRingInd;
VT_img = objectiveStruct.VT_img;

appliedPressureSys = objectiveStruct.appliedPressureSys ;
appliedPressureDias = objectiveStruct.appliedPressureDias ;

febio_spec=objectiveStruct.febio_spec ;
febioAnalysis= objectiveStruct.febioAnalysis;
febioFebFileName= objectiveStruct.febioFebFileName;
savePath= objectiveStruct.savePath;
febioLogFileName_disp = objectiveStruct.febioLogFileName_disp;
febioLogFileName_stress=objectiveStruct.febioLogFileName_stress;
U_find= objectiveStruct.U;


%% Get optimization parameters
K=S;
C=Pn;

%% Creating undeformed state
V_unDEF = VT_img - K*U_find;

%%  Inflate till diastole
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_unDEF,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V_unDEF; %The nodel coordinates
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V_unDEF,1);

febio_spec.Material.material{1}.c1= C; %setting new stiffness
febio_spec.Material.material{1}.c2=C;

%Inflate till diastole
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureDias];

%Exporting the FEBio input file
febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%START FEBio
[runFlag]=runMonitorFEBio(febioAnalysis);


if runFlag==1
    
    %Importing displacement
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    N_disp_mat=dataStruct.data; %Displacement
    U_dias = N_disp_mat(:,:,end);
    
    
    %Derive deformed configuration
    V_DEF_dias= V_unDEF + U_dias;
    
    %Obtained systolic outer ring and mean radius
    V_DEF_dias_innerRing = V_DEF_dias(indVer,:);
    mean_V_DEF_dias = mean(V_DEF_dias_innerRing);
    allRadialVecs_dias = V_DEF_dias_innerRing - (mean_V_DEF_dias).*ones (size(V_DEF_dias_innerRing));
    meanRadius_dias = mean(sqrt(sum( allRadialVecs_dias.^2,2)));
    
end


%% Inflate till systole
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_unDEF,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V_unDEF; %The nodel coordinates
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V_unDEF,1);

febio_spec.Material.material{1}.c1= C; %setting new stiffness
febio_spec.Material.material{1}.c2=C;
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureSys];

%Exporting the FEBio input file
febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%START FEBio
[runFlag]=runMonitorFEBio(febioAnalysis);


if runFlag==1
    
    %Importing displacement
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    N_disp_mat=dataStruct.data; %Displacement
    U_sys = N_disp_mat(:,:,end);
    %import elemental stresses from logfile
    dataStruct_stress = importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
    E_stress_mat= dataStruct_stress.data; %Elemental stresses   

    
    %Derive deformed configuration
    V_DEF_sys= V_unDEF + U_sys;
    %Obtained systolic outer ring and mean radius
    V_DEF_sys_innerRing = V_DEF_sys(indVer,:);
    mean_V_DEF_sys = mean(V_DEF_sys_innerRing);
    allRadialVecs_sys = V_DEF_sys_innerRing - (mean_V_DEF_sys).*ones (size(V_DEF_sys_innerRing));
    meanRadius_sys = mean(sqrt(sum( allRadialVecs_sys.^2,2)));
    
    Popt = abs((r_exp_Sys - meanRadius_sys)-(r_exp_Dias-meanRadius_dias));
    
    OPT_stats_out.Popt=Popt;
    
    
else
    Popt=NaN(2,1); %Squared differences
    OPT_stats_out=[];
end
end
