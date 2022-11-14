%% ckOptimize_aorta2stepOptimize_rootDisp

tic
clear; close all; clc;

%% Control parameters
smoothFactorCentreLine=1;
pointSpacingLoft = 8;
numSplitSteps=2;

bezierTangency=0.5;
numSmoothBranchAttachments=50;
smoothIncludeSteps=4;


f=1/3;
smoothParLoftCurve=0.2;
numThickenSteps =1;
nodalThickness= 2;

optionStructRayTrace.tolEps        = 1e-6;
optionStructRayTrace.triSide       = 1;
optionStructRayTrace.rayType       = 'ray';
optionStructRayTrace.exclusionType = 'inclusive';
optionStructRayTrace.paired        = 0;

%% Plot settings

fontSize=15;

%Define plot color (black or white background)
colorMode=1;
switch colorMode
    case 1
        figStruct.ColorDef='white';
        figStruct.Color='w';
    case 2
        figStruct.ColorDef='black';
        figStruct.Color='k';
end

%% FE control parameters

% Path names
savePath = fileparts(mfilename('fullpath'));

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt'];%Log file name for exporting stress

%Load
hG2N = 0.000133322;
pressureDiasmmHg = 80;
pressureSysmmHg = 130;

appliedPressureDias= -(pressureDiasmmHg*hG2N);
appliedPressureSys= -(pressureSysmmHg*hG2N);



displacementMagnitude=8;
%Longitudinal diplacement for inflation till diastolic phase systolic phase
dispDias = [0 0 0];
dispSys = [ 5.7   -1.7   -6.4]/appliedPressureSys; 

%Material parameter set
c1= 1.1; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
% c1=0.00005; %Shear-modulus-like parameter
% m1=1; %Material parameter setting degree of non-linearity
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

%Initial tol for parameter optimization
tol = 10;

%% importlumen and centreline path
%Load MRI data
load('MRI_extracted');

%% Main trunk
% Lumen
output_lumen = V_MRI_extract;

output_lumen{1, 1} = flip(output_lumen{1, 1});
output_lumen{1, 5} = flip(output_lumen{1, 5});
output_lumen{1, 6} = flip(output_lumen{1, 6});
%Center line path
output_centre = Centre;

%% Branch data

numPointsEllipse=100;

V_BCA_Br1_Ori =xlsread('BM1_level1');
V_BCA_Br1_Ori_e=fitEllipse(V_BCA_Br1_Ori,numPointsEllipse);

V_BCA_Br1=xlsread('BM1_level2');
V_BCA_Br1_e=fitEllipse(V_BCA_Br1,numPointsEllipse);

V_BCA_Br2_Ori =  (xlsread('BM2_level1'));
V_BCA_Br2_Ori_e=fitEllipse(V_BCA_Br2_Ori,numPointsEllipse);

V_BCA_Br2=  ((xlsread('BM2_level2')));
V_BCA_Br2_e=fitEllipse(V_BCA_Br2,numPointsEllipse);

V_BCA_Br3_Ori = xlsread('BM3_level1');
V_BCA_Br3_Ori_e=fitEllipse(V_BCA_Br3_Ori,numPointsEllipse);

V_BCA_Br3=xlsread('BM3_level2');
V_BCA_Br3_e=fitEllipse(V_BCA_Br3,numPointsEllipse);


%%
cFigure(figStruct);hold on;

plotV(V_BCA_Br1_Ori,'r.','markerSize',25);
plotV(V_BCA_Br1_Ori_e,'k-','LineWidth',2);

plotV(V_BCA_Br2_Ori,'g.-','markerSize',25);
plotV(V_BCA_Br2_Ori_e,'k-','LineWidth',2);

plotV(V_BCA_Br3_Ori,'b.-','markerSize',25);
plotV(V_BCA_Br3_Ori_e,'k-','LineWidth',2);

plotV(V_BCA_Br1,'r.-','markerSize',25);
plotV(V_BCA_Br1_e,'k-','LineWidth',2);

plotV(V_BCA_Br2,'g.-','markerSize',25);
plotV(V_BCA_Br2_e,'k-','LineWidth',2);

plotV(V_BCA_Br3,'b.-','markerSize',25);
plotV(V_BCA_Br3_e,'k-','LineWidth',2);

axisGeom(gca,fontSize); camlight headlight;
drawnow;

%% Extract centreline from raw data (main trunk + lumen)
V_cent = output_centre;
V_cent_BCA_Br1= [mean(V_BCA_Br1_Ori);mean(V_BCA_Br1)];
V_cent_BCA_Br2= [mean(V_BCA_Br2_Ori);mean(V_BCA_Br2)];
V_cent_BCA_Br3= [mean(V_BCA_Br3_Ori);mean(V_BCA_Br3)];

%% Resampling aorta section contours and smoothing
d=zeros(size(output_lumen,2),1);

for indNow=1:1:size(output_lumen,2)
    V_lumen= cell2mat(output_lumen{1,indNow});
    d(indNow)=max(pathLength([V_lumen;V_lumen(1,:)]));
end
nSegment=round(max(d)/pointSpacingLoft);


contourSmooth = output_lumen;
for indNow=1:1:size(output_lumen,2)
    %Resample section contour
    V_lumen_original=cell2mat(output_lumen{1,indNow}); %Current contour
    d=pathLength(V_lumen_original);
    V_lumen_smooth=evenlySampleCurve(V_lumen_original,nSegment,0.01,1); %Resample evenly
    contourSmooth{1,indNow} = V_lumen_smooth;
end

%% Offsetting section curves outward if thickening is inward
contourSmooth_pre = contourSmooth;

for q=1:size(contourSmooth,2)
    Vc=contourSmooth{1,q}; %Current curve vertices
    contourSmooth{1,q}=curveOffset(Vc,nodalThickness);
end


%% Centreline of aorta Bezier curve/cubic polynomial

for q=1:1:size(contourSmooth,2)
    
    p1.centre{q}= mean(contourSmooth_pre{1,q},1); %Curve center
    vf=-vecnormalize(contourSmooth_pre{1,q}-[contourSmooth_pre{1,q}(2:end,:);contourSmooth_pre{1,q}(1,:)]); %Allong curve path vectors
    vb=vecnormalize(contourSmooth_pre{1,q}-[contourSmooth_pre{1,q}(end,:);contourSmooth_pre{1,q}(1:end-1,:)]); %Allong curve path vectors
    v=(vf+vb)/2;
    r=vecnormalize(contourSmooth_pre{1,q}-p1.centre{q}); %Position vector wrt mean
    v1=vecnormalize(cross(v,r)); %perimeter quasi Z-vectors
    n.normal{q} = vecnormalize(mean(v1,1));
end

cFigure(figStruct); hold on;
for q=1:1:(size(contourSmooth,2)-1)
    
    v1 = contourSmooth{1,q};
    v2 = contourSmooth{1,q+1};
    
    cp1 =  p1.centre{q};
    cp2 =  p1.centre{q+1};
    n1= n.normal{q};
    n2= n.normal{q+1};
    
    
    numPoints = 100;
    Vg=sweepCurveSmooth(cp1,cp2,n1,n2,numPoints,smoothParLoftCurve,f);
%         Vg= evenlySampleCurve(Vg,20,'pchip',1);
    Vcentre.centre{q} = Vg;
    
    % View plots
    hp(1)=plotV(cp1,'r.','MarkerSize',50);
    hp(2)=quiverVec(cp1,n1,10,'r');
    hp(3)=plotV(cp2,'g.','MarkerSize',50);
    hp(4)=quiverVec(cp2,n2,10,'g');
    hp(5)=plotV(v1,'r.-','LineWidth',3,'MarkerSize',15);
    hp(6) = plotV(v2,'g.-','LineWidth',3,'MarkerSize',15);
    hp(5) = plotV(Vcentre.centre{q},'k.-','LineWidth',3,'MarkerSize',25);
    axisGeom(gca,fontSize);
    
end


%% Viewing lumen-lumen connecting lines on original contours
cFigure(figStruct);hold on;
for nLum = 1:1:(size(contourSmooth,2)-1)
    V1 = contourSmooth{1,nLum};
    V2 = contourSmooth{1,nLum+1};
    d = getInterPoints(V1,V2);
    
    %Before refinement plot of one segment
    
    for i = 1:1:size(V1,1)
        plotV(V1,'r.-','Linewidth',3)
        plotV(V2,'g.-','Linewidth',3)
        plotV(d.p{1,i},'k.-','Markersize',25)
        hold on;
    end
    axisGeom(gca,fontSize);
    drawnow;
    
end

%% Minimizing lumen-lumen connecting lines distances
V2_find = contourSmooth{1,1}; %define V2_find with first original fixed lumen
for nLum = 1:1:(size(contourSmooth,2)-1)
    
    V1 = V2_find ;
    V2 = contourSmooth{1,nLum+1};
    d = getInterPoints(V1,V2);
    
    %Minimizing inter point distances
    indDist = zeros(size(V1,1),1); % distance between two points on two different lumen
    sumDist = zeros(size(V1,1),1); %sum of distances between all points
    
    
    for q = 1:1:size(V1,1)
        
        for i = 1:1:size(V1,1)
            a = d.p{1,i}(1,:,:);
            b = d.p{1,i}(2,:,:);
            indDist(i) = norm(a-b);
        end
        
        V_find.index{q} =V2; %store all the refined lumens in V_find
        sumDist(q) = sum(indDist);
        V2 = [V2(2:end,:,:);V2(1,:,:)];
        d = getInterPoints(V1,V2);
    end
    
    
    index = find(sumDist==min(sumDist));
    V2_find = V_find.index{index};
    contourSmoothUpdate{nLum+1}= V_find.index{index}; %storing updated contours
end


%% Viewing connecting lines after minimizing distance
contourSmoothUpdate{1,1}= contourSmooth{1,1}; %filling up the first empty cell of updated structure
cFigure(figStruct); hold on;
for nLum = 1:1:(size(contourSmoothUpdate,2)-1)
    
    V1 = contourSmoothUpdate{1,nLum};
    V2 = contourSmoothUpdate{1,nLum+1};
    d = getInterPoints(V1,V2);
    
    for i = 1:1:size(V1,1)
        plotV(V1,'r.-','Linewidth',3)
        plotV(V2,'g.-','Linewidth',3)
        plotV(d.p{1,i},'k.-','Markersize',25) %plotting lines between two lumen
        hold on;
    end
    axisGeom(gca,fontSize);
    drawnow;
end


%% Lofting

cFigure(figStruct); hold on;
axisGeom(gca,fontSize);
camlight headlight;

maxC=0;
numLofts=(size(contourSmoothUpdate,2)-1);
F_main_cell=cell(numLofts,1);
V_main_cell=cell(numLofts,1);
C_main_gradient_cell=cell(numLofts,1);
C_main_feature_cell=cell(numLofts,1);
indStart_cell=cell(numLofts,1);
indEnd_cell=cell(numLofts,1);
for q=1:1:numLofts
    
    V1 = contourSmoothUpdate{1,q};
    V2 = contourSmoothUpdate{1,q+1};
    n1 = n.normal{q};
    n2 = n.normal{q+1};
    Vg= Vcentre.centre{1,q};
    d= mean([max(pathLength([V1;V1(1,:)])) max(pathLength([V2;V2(1,:)]))]);
    pointSpacingLong = d/nSegment;
    
    Vg = evenlySpaceCurve (Vg, pointSpacingLong, 'pchip');
    %     Vg = evenlySampleCurve (Vg, pointSpacingLong, 'pchip');
    Vcentre.centre{1,q} =Vg;
    numPointsLong = size(Vg,1);
    
    
    Xc = zeros(numPointsLong,size(V1,1));
    Yc = zeros(numPointsLong,size(V1,1));
    Zc = zeros(numPointsLong,size(V1,1));
    f_seg=f;
    for qc = 1:1:size(V1,1)
        %         Vc=sweepCurveBezier(V1(qc,:),V2(qc,:),n1,n2,numPointsLong,f_seg); %Compute curve
        Vc=sweepCurveSmooth(V1(qc,:),V2(qc,:),n1,n2,numPointsLong,f); %Compute curve
        Xc(:,qc)=Vc(:,1);
        Yc(:,qc)=Vc(:,2);
        Zc(:,qc)=Vc(:,3);
    end
    %Color data
    c=(1:1:numPointsLong)';
    C=c(:,ones(1,size(V1,1))); %Color gradient for each vertex (no. of points along length of curve x no. of points along circumference of the curve)
    
    %Create quad patch data
    [F,V,C] = surf2patch(Xc,Yc,Zc,C);% C obtained from this step is still applied to vertices
    
    %Close section if required
    
    I=[(1:size(Zc,1)-1)' (1:size(Zc,1)-1)' (2:size(Zc,1))' (2:size(Zc,1))' ];
    J=[size(Zc,2).*ones(size(Zc,1)-1,1) ones(size(Zc,1)-1,1) ones(size(Zc,1)-1,1) size(Zc,2).*ones(size(Zc,1)-1,1)];
    F_sub=sub2ind(size(Zc),I,J);
    F=[F;F_sub];
    C(end-size(F_sub,1):end,:)=C(end-size(F_sub,1):end,:)+0.5;
    [C]=vertexToFaceMeasure(F,C); %Convert vertex colors to face colors
    C=round(C)-1; % here size(C,1)= no. of faces
    
    %Get indices for start and end curves aka rings
    indStart=1:numPointsLong:size(V,1);
    indEnd=numPointsLong:numPointsLong:size(V,1);
    
    F_main_cell{q}=F;
    V_main_cell{q}=V;
    C_main_gradient_cell{q} = (C)+ maxC;
    C_main_feature_cell{q} = q*ones(size(F,1),1);
    maxC=maxC+max(C);
    indStart_cell{q}=indStart;
    indEnd_cell{q}=indEnd;
    
    % Plotting geometry
    h(1) = plotV(Vg,'k.-','LineWidth',3);
    h(2) = plotV(V(indStart,:), 'r.-','Markersize',25,'LineWidth',3);
    h(3) = plotV(V(indEnd,:), 'g.-','Markersize',35,'LineWidth',2);
    
    h=gpatch(F,V,C,'k',0.85);
    legend(h,{'Lofted surface'});
    drawnow;
    
end

%% Merging lofted features

[F_main,V_main,C_main]=joinElementSets(F_main_cell,V_main_cell,C_main_gradient_cell);
[~,~,C_main_feature]=joinElementSets(F_main_cell,V_main_cell,C_main_feature_cell);
F_main=fliplr(F_main); %Invert orientation

[F_main,V_main,ind1,indFix]=mergeVertices(F_main,V_main);
Eb = patchBoundary(F_main);

%Correcting ring indices
indexOffset=0;
for q=1:1:numLofts
    
    %Correct for joining if nodes for 2nd and further sections
    if q>1
        indexOffset=indexOffset+size(V_main_cell{q-1},1); %Add current number of nodes to the offset
        indStart_cell{q}=indStart_cell{q}+indexOffset; %Correct for the joining of nodes
        indEnd_cell{q}=indEnd_cell{q}+indexOffset; %Correct for the joining of nodes
    end
    
    %Correct for the merging nodes
    indStart_cell{q}=indFix(indStart_cell{q});
    indEnd_cell{q}=indFix(indEnd_cell{q});
end

%%

cFigure(figStruct); hold on;
gpatch(F_main,V_main,C_main_feature,'k',0.85);

for q=1:1:numLofts
    plotV(V_main(indStart_cell{q},:), 'r.-','Markersize',25,'LineWidth',3);
    plotV(V_main(indEnd_cell{q},:), 'g.-','Markersize',35,'LineWidth',2);
end

axisGeom(gca,fontSize);
camlight headlight;
colormap gjet; icolorbar;
drawnow;

%% Mesh refinement

% [F_main,V_main,C_sub,~]=subQuad(F_main,V_main,numSplitSteps,1);
[F_main,V_main,C_sub]=subQuadCatmullClark(F_main,V_main,numSplitSteps,0);

C_main = C_main(C_sub); %expanding color for the refined mesh
C_main_feature = C_main_feature(C_sub); %expanding color for the refined mesh

%%
cFigure(figStruct);hold on;
gpatch(F_main,V_main,C_main_feature,'k',0.85)
axisGeom(gca,fontSize);
camlight headlight;
colormap gjet; icolorbar;
gdrawnow;

%% Add branches

cutDist=mean(patchEdgeLengths(F_main,V_main));

% C_main_feature=zeros(size(C_main));
[F_main,V_main,C_main,C_main_feature]=cutMesh(F_main,V_main,C_main,C_main_feature,V_BCA_Br1_e,V_BCA_Br1_Ori_e,optionStructRayTrace,cutDist,bezierTangency,numSmoothBranchAttachments,smoothIncludeSteps);
[F_main,V_main,C_main,C_main_feature]=cutMesh(F_main,V_main,C_main,C_main_feature,V_BCA_Br2_e,V_BCA_Br2_Ori_e,optionStructRayTrace,cutDist,bezierTangency,numSmoothBranchAttachments,smoothIncludeSteps);
[F_main,V_main,C_main,C_main_feature]=cutMesh(F_main,V_main,C_main,C_main_feature,V_BCA_Br3_e,V_BCA_Br3_Ori_e,optionStructRayTrace,cutDist,bezierTangency,numSmoothBranchAttachments,smoothIncludeSteps);

%%

cFigure(figStruct);hold on;
% gpatch(F_main,V_main,C_main,'k',1)

gpatch(F_main,V_main,C_main,'k',1)
% patchNormPlot(F_main,V_main);

plotV(V_BCA_Br1_Ori_e,'r.-','markerSize',10);
plotV(V_BCA_Br1_e,'r.-','markerSize',10);

plotV(V_BCA_Br2_Ori_e,'g.-','markerSize',10);
plotV(V_BCA_Br2_e,'g.-','markerSize',10);

plotV(V_BCA_Br3_Ori_e,'b.-','markerSize',10);
plotV(V_BCA_Br3_e,'b.-','markerSize',10);

axisGeom(gca,fontSize); camlight headlight;
colormap gjet; colorbar;
drawnow;

%%
% plot

cFigure(figStruct); hold on;
gpatch(F_main,V_main,C_main_feature,'k',1);
% for q=1:1:numel(V_endCurve_cell)
%     plotV(V_endCurve_cell{q},'m.-','markerSize',25);
% end
axisGeom(gca,fontSize);
colormap gjet; icolorbar;
camlight headlight;
drawnow;

%% Thicken to create hexahedral elements

[ET,VT,Fp1,Fp2]=patchThick(F_main,V_main,-1,nodalThickness,numThickenSteps);
ET=ET(:,[5:8 1:4]); %Invert elements due to thicken direction
C_ET_main=repmat(C_main,numThickenSteps,1);
C_ET_feature=repmat(C_main_feature,numThickenSteps,1);
FT_inner=Fp2;
FT_outer=Fp1;
FT_inner = fliplr(FT_inner);

%Get element volumes
[vol_E,logicPositive]=hexVol(ET,VT); %Absolute volumes
vol_E(~logicPositive)=-vol_E(~logicPositive); %Flip negative volumes to show up as negative

[FT,mapElementFaceData,C_FT_type]=element2patch(ET);
C_FT_main=C_ET_main(mapElementFaceData);
C_FT_feature=C_ET_feature(mapElementFaceData);
C_FT_vol_E=vol_E(mapElementFaceData);


%Get boundary faces
indBoundaryFaces = tesBoundary(FT);
FT_b= FT(indBoundaryFaces,:);
C_FT_b_main=C_FT_main(indBoundaryFaces,:);
C_FT_b_type=C_FT_type(indBoundaryFaces,:);
C_FT_b_feature=C_FT_feature(indBoundaryFaces);
C_FT_b_vol_E=C_FT_vol_E(indBoundaryFaces);

%%

cFigure(figStruct);
subplot(1,2,1); hold on;
gpatch(FT_b,VT,C_FT_b_main,'k',1);
axisGeom(gca,fontSize); camlight headlight;
colormap gjet; colorbar;

subplot(1,2,2);hold on;
gpatch(FT_b,VT,C_FT_b_feature,'k',1);
axisGeom(gca,fontSize); camlight headlight;
colormap gjet; icolorbar;
drawnow



%% Get indices for nodes that are a member of the branchces
logicFacesExclude=C_FT_b_feature>numLofts;
indexVerticesExclude=unique(FT_b(logicFacesExclude,:));

%Inner nodes
indInnerRingNodes=[];
indOuterRingNodes=[];
for q=1:1:numLofts
    
    %Inner faces/edge/vertices
    logicNow = (C_FT_b_type==1) & (C_FT_b_feature==q); %Logic for current faces
    Fb_now=FT_b(logicNow,:); %Current face select
    Eb_Now=patchBoundary(Fb_now); %All boundary edges of the current face set
    logicKeep=all(~ismember(Eb_Now,indexVerticesExclude),2);
    Eb_Now=Eb_Now(logicKeep,:);
    indInnerRingNodes=[indInnerRingNodes; Eb_Now(:)];
    indInnerRingNodes=unique(indInnerRingNodes);
    
    %Outer faces/edge/vertices
    logicNow = (C_FT_b_type==2) & (C_FT_b_feature==q); %Logic for current faces
    Fb_now=FT_b(logicNow,:); %Current face select
    Eb_Now=patchBoundary(Fb_now); %All boundary edges of the current face set
    logicKeep=all(~ismember(Eb_Now,indexVerticesExclude),2);
    Eb_Now=Eb_Now(logicKeep,:);
    indOuterRingNodes=[indOuterRingNodes; Eb_Now(:)];
    indOuterRingNodes=unique(indOuterRingNodes);
end

%% Get inner ring nodes at ascending aorta level

logicFacesSegment1=C_FT_b_feature==1;
logicFacesSegment2=C_FT_b_feature==2;
indexVertices_1=unique(FT_b(logicFacesSegment1,:));
indexVertices_2=unique(FT_b(logicFacesSegment2,:));

logicInner_12=ismember(indInnerRingNodes,indexVertices_1) & ismember(indInnerRingNodes,indexVertices_2);
logicOuter_12=ismember(indOuterRingNodes,indexVertices_1) & ismember(indOuterRingNodes,indexVertices_2);

tNow=sqrt(sum((VT(indInnerRingNodes(logicInner_12),:)-VT(indOuterRingNodes(logicOuter_12),:)).^2,2));


indRingAsc= indInnerRingNodes(logicInner_12);
%%
cFigure(figStruct); hold on;
gpatch(FT_b,VT,C_FT_b_feature,'none',0.25);
% gpatch(FT_b(logicNow,:),VT,'kw','k',1);
% gpatch(Eb_Now,VT,'none','r',1,3);
plotV(VT(indRingAsc,:),'r.-','MarkerSize',25);
plotV(VT(indOuterRingNodes(logicOuter_12),:),'b.-','MarkerSize',25);
axisGeom(gca,fontSize); camlight headlight;
colormap(gca,gjet(250)); icolorbar([],gca);
gdrawnow;


%%

nBins=25;
[nHist,eHist] = histcounts(vol_E,nBins);
w=mean(diff(eHist));
xHist=eHist(1:end-1)+w/2;

maxQ=max(abs(vol_E));

cFigure(figStruct);
gtitle('Element volumes',fontSize)
subplot(1,2,1);hold on;
gpatch(FT_b,VT,C_FT_b_vol_E,'k',1);
axisGeom(gca,fontSize); camlight headlight;
colormap spectral; colorbar; caxis([-maxQ maxQ]);

subplot(1,2,2);hold on;
h=bar(xHist,nHist,w*2,'FaceColor','flat');
h.CData = xHist;
colormap spectral; caxis([-maxQ maxQ]);
set(gca,'FontSize',fontSize);
drawnow

%%

if any(logicPositive==0)
    warning('Negative hex volume found');
end

%%
cFigure(figStruct);
% subplot(1,3,1); hold on;
% title('Face types');
gpatch(FT_b,VT,C_FT_b_type,'k');
axisGeom(gca,fontSize); camlight headlight;
colormap(gca,gjet(250)); %icolorbar([],gca);

subplot(1,3,2); hold on;
title('Main trunk path');
gpatch(FT_b,VT,C_FT_b_main,'k');
axisGeom(gca,fontSize); camlight headlight;
colormap(gca,gjet(250)); colorbar(gca);

subplot(1,3,3); hold on;
title('Feature label');
gpatch(FT_b,VT,C_FT_b_feature,'k');
axisGeom(gca,fontSize); camlight headlight;
colormap(gca,gjet(250)); icolorbar([],gca);

gdrawnow;

%% Define boundary conditions

bcSupportList=unique(FT_b( ismember(C_FT_b_type,[3 5]),:));
bcPrescribeList=unique(FT_b(ismember(C_FT_b_type,4),:));

F_pressure=fliplr(FT_inner);

N_disp=mean(patchNormal(FT_b(ismember(C_FT_b_type,4),:),VT),1);
U_applied = displacementMagnitude * N_disp;
%%

cFigure(figStruct); hold on;
gpatch(FT_b,VT,'w','none',0.5);

hp1=gpatch(F_pressure,VT,'bw','k');
patchNormPlot(F_pressure,VT);
hp2=plotV(VT(bcSupportList,:),'k.','MarkerSize',25);
hp3=plotV(VT(bcPrescribeList,:),'r.','MarkerSize',25);
% hp4=plotV(VT(indInnerRingNodes(logicInner_12),:),'g.-','MarkerSize',25);

hp4=quiverVec(VT(bcPrescribeList,:),U_applied(ones(numel(bcPrescribeList),1),:),[],'r');

legend([hp1 hp2 hp3],{'Pressure surface','Supported nodes','Prescribed displacement nodes','Displacement vectors'})
% legend([hp1 hp2 hp3],{'Pressure surface','Supported nodes'})
axisGeom(gca,fontSize);
camlight headlight;
gdrawnow;


VT_img =VT; %store vertices of orginal geometry


%% Acquiring diastolic and systolic nodal displacmenet wrt undeformed state

appliedPressure = appliedPressureSys;


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
% febio_spec.Mesh.Elements{1}.elem.scale.ATTR.lc=1; %The element matrix
% febio_spec.Mesh.Elements{1}.elem.scale.VAL=Estress; %The element matrix

% -> Surfaces
surfaceName1='LoadedSurface';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure;

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

nodeSetName2='bcPrescribeList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);


%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.type='fix';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dofs='x,y,z';

% displacementMagnitude = [ 5.7   -1.7   -6.4];
displacementMagnitude = dispDias;
febio_spec.Boundary.bc{2}.ATTR.type='prescribe';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{2}.dof='x';
febio_spec.Boundary.bc{2}.scale.ATTR.lc=1;
febio_spec.Boundary.bc{2}.scale.VAL=displacementMagnitude(1);
febio_spec.Boundary.bc{2}.relative=0;

febio_spec.Boundary.bc{3}.ATTR.type='prescribe';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{3}.dof='y';
febio_spec.Boundary.bc{3}.scale.ATTR.lc=1;
febio_spec.Boundary.bc{3}.scale.VAL=displacementMagnitude(2);
febio_spec.Boundary.bc{3}.relative=0;

febio_spec.Boundary.bc{4}.ATTR.type='prescribe';
febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{4}.dof='z';
febio_spec.Boundary.bc{4}.scale.ATTR.lc=1;
febio_spec.Boundary.bc{4}.scale.VAL=displacementMagnitude(3);
febio_spec.Boundary.bc{4}.relative=0;

%Loads section
% -> Surface load
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
% febio_spec.Loads.surface_load{1}.pressure.VAL= appliedPressure;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=1;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 appliedPressure];

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
febioAnalysis.maxLogCheckTime=60000; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!




%% Import FEBio results

if runFlag==1 %i.e. a succesful run
    
    %Importing displacement
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    N_disp_mat=dataStruct.data; %Displacement
    U = N_disp_mat(:,:,end);  
    
    
    %import elemental stresses from logfile
    dataStruct_stress = importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
    E_stress_mat= dataStruct_stress.data; %Elemental stresses
end


U_find = U;

VT_MRI = VT_img;

%% Creating experimental undeformed state and finding synthetic diastolic and sytolic configurations

V_unDEF = VT_MRI - 0.3*U_find;

indVer=indRingAsc;
V_unDEF_innerRing = V_unDEF(indVer,:);
mean_V_unDEF = mean(V_unDEF_innerRing);
allRadialVecs_unDEF = V_unDEF_innerRing - (mean_V_unDEF).*ones (size(V_unDEF_innerRing));
meanRadius_unDEF = mean(sqrt(sum( allRadialVecs_unDEF.^2,2)));

%% Run simulations in two steps
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_unDEF,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V_unDEF; %The nodel coordinates
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V_unDEF,1);

% febio_spec.Material.material{1}.c1= c1; %setting new stiffness

%Inflate till diastole
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureDias];

%Boundary conditions for diastole
febio_spec.Boundary.bc{2}.scale.VAL=dispDias(1);
febio_spec.Boundary.bc{3}.scale.VAL=dispDias(2);
febio_spec.Boundary.bc{4}.scale.VAL=dispDias(3);

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
    
    %Derive deformed configuration
    V_DEF_dias= V_unDEF + U_dias;
    VT_img = V_DEF_dias; %redefining VT_img
    
    %Obtained systolic outer ring and mean radius
    V_DEF_dias_innerRing = V_DEF_dias(indVer,:);
    mean_V_DEF_dias = mean(V_DEF_dias_innerRing);
    allRadialVecs_dias = V_DEF_dias_innerRing - (mean_V_DEF_dias).*ones (size(V_DEF_dias_innerRing));
    r_exp_Dias = mean(sqrt(sum( allRadialVecs_dias.^2,2))); %experimental diastolic state
    
end


%Inflate till systole
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureSys];

%Boundary conditions for systole
febio_spec.Boundary.bc{2}.scale.VAL=dispSys(1);
febio_spec.Boundary.bc{3}.scale.VAL=dispSys(2);
febio_spec.Boundary.bc{4}.scale.VAL=dispSys(3);


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
    r_exp_Sys = mean(sqrt(sum( allRadialVecs_sys.^2,2)));
end

% cFigure; hold on;
% gpatch(FT,VT,'rw','none',1);
% gpatch(FT,V_DEF_sys,'none','k',1);
% plotV(V_DEF_sys,'g.','Markersize',20)
% plotV(V_DEF_sys(bcPrescribeList,:),'r.','Markersize',20)
% quiverVec(V_DEF_sys(bcPrescribeList,:),U_applied(ones(numel(bcPrescribeList),1),:),[],'r');
% axis off
% axisGeom





%% 
cFigure;
hold on;
H_opt = plot(0,NaN,'r.-','Markersize',25);
view(2); axis tight;  grid on;
set(gca,'FontSize',fontSize);
xlabel('Iteration no.');
ylabel('C val');

cFigure;
hold on;
O_opt = plot(0,NaN,'g.-','Markersize',25);
view(2); axis tight;  grid on;
set(gca,'FontSize',fontSize);
xlabel('Iteration no.');
ylabel('Objective function val');

% cFigure;
% hold on;
% H_opt(2) = plot(0,NaN,'b.-','Markersize',25);
% view(2); axis tight;  grid on;
% set(gca,'FontSize',fontSize);
% xlabel('Iteration no.');
% ylabel('Parameter value');


drawnow;

%% Creating Objetcive struct for getting optimization and determining U_find
objectiveStruct.hOpt=H_opt;
objectiveStruct.oOpt=O_opt;
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
%0.8278


%% Start parameter loop here
% counter = 1;
% kinit = 1;
% Cinit= 0.5;
% res= 1;
% res= [];
K_new = [];

% cFigure; hold on;
%% Start parameter loop here
counter = 1;
kinit = 0.25;
Cinit= 1.5;
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
    
    
    
    %% Get material parameter for a given shrinkage factor
    
    Pn = C;    
    [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveBackwardInc(Pn,objectiveStruct,r_exp_Sys,r_exp_Dias,S,VT_MRI,dispDias,dispSys),Pn,[],[],OPT_options);    

    C = Pn_opt;

    %% Get shrinkage factor for a given material parameter
    Sn= S;
    [Sn_opt,SPT_out.resnorm,SPT_out.residual]= lsqnonlin(@(Sn) getShrinkFact(Sn,objectiveStruct,C,dispDias,VT_MRI),Sn,[],[],SPT_options);    
  
    S= Sn_opt;
    
    %% Creating undeformed state
    V_unDEF = VT_MRI - S*U_find;
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
    febio_spec.Material.material{1}.c2= C; %setting new stiffness
    
    %Inflate till diastole
    febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureDias];
    
    %Boundary conditions for diastole
    febio_spec.Boundary.bc{2}.scale.VAL=dispDias(1);
    febio_spec.Boundary.bc{3}.scale.VAL=dispDias(2);
    febio_spec.Boundary.bc{4}.scale.VAL=dispDias(3);
    
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
    
    %Boundary conditions for diastole
    febio_spec.Boundary.bc{2}.scale.VAL=dispSys(1);
    febio_spec.Boundary.bc{3}.scale.VAL=dispSys(2);
    febio_spec.Boundary.bc{4}.scale.VAL=dispSys(3);
    
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
    
    %     res = abs((r_exp_Sys - meanRadius_sys ) + (r_exp_Dias-meanRadius_dias));
    res = abs((r_exp_Sys - meanRadius_sys ) - (r_exp_Dias - meanRadius_dias));
    counter= counter+1;
end


save aortaOptimize8.mat  

toc


%% FUNCTIONS
function [Sopt,SPT_stats_out]=getShrinkFact(Sn,objectiveStruct,C,dispDias,VT_MRI)

%% Get variables
% H_opt = objectiveStruct.hOpt;
VT_img = objectiveStruct.VT_img;
appliedPressureDias = objectiveStruct.appliedPressureDias ;
febio_spec=objectiveStruct.febio_spec;
febioAnalysis= objectiveStruct.febioAnalysis;
febioFebFileName= objectiveStruct.febioFebFileName;
savePath= objectiveStruct.savePath;
febioLogFileName_disp = objectiveStruct.febioLogFileName_disp;
indRingAsc=objectiveStruct.innerRingInd;
r_exp_Dias=objectiveStruct.r_exp_Dias  ;

% P= Pn.*objectiveStruct.parNormFactors;
U_find= objectiveStruct.U;


%% Get optimization parameters
K = Sn;

%% Creating undeformed state
V_unDEF = VT_MRI - K*U_find;


%% Inflate till diastole
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_unDEF,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V_unDEF; %The nodel coordinates
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V_unDEF,1);

febio_spec.Material.material{1}.c1= C; %setting new stiffness
febio_spec.Material.material{1}.c2= C; %setting new stiffness

febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureDias];

%Boundary conditions for diastole
%Boundary conditions for diastole
febio_spec.Boundary.bc{2}.scale.VAL=dispDias(1);
febio_spec.Boundary.bc{3}.scale.VAL=dispDias(2);
febio_spec.Boundary.bc{4}.scale.VAL=dispDias(3);

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
    
    

    indVer=indRingAsc;
    %Obtained systolic outer ring and mean radius
    V_DEF_dias_innerRing = V_DEF_dias(indVer,:);
    mean_V_DEF_dias = mean(V_DEF_dias_innerRing);
    allRadialVecs_dias = V_DEF_dias_innerRing - (mean_V_DEF_dias).*ones (size(V_DEF_dias_innerRing));
    meanRadius_dias = mean(sqrt(sum( allRadialVecs_dias.^2,2)));



    V_unDEF_ir = V_unDEF(indVer,:);
    mean_V_unDEF_ir = mean( V_unDEF_ir);
    allRadialVecs_ir =  V_unDEF_ir - ( mean_V_unDEF_ir).*ones (size(V_unDEF_ir));
    meanRadius_ir = mean(sqrt(sum(  allRadialVecs_ir.^2,2)));

    vDev = VT_img - V_DEF_dias;
    Sopt = (sqrt(sum(vDev.^2,2)));
    %     Sopt = vDev;

    %     Sopt = abs(r_exp_Dias - meanRadius_dias);
    
    SPT_stats_out.Sopt=Sopt;
%     if ~isempty(H_opt)
% %            H_opt(1).XData= [H_opt(1).XData H_opt.XData(end)+1];
% %            H_opt(1).YData= [H_opt(1).YData (Popt)]; 
%            H_opt.XData= [H_opt.XData H_opt.XData(end)+1];
%            H_opt.YData= [H_opt.YData Sn]; 
%            
%            drawnow;
%        end
%     
    
    
    
else
    Sopt=NaN(2,1); %Squared differences
    SPT_stats_out=[];
end
end


function [Popt,OPT_stats_out]=objectiveBackwardInc(Pn,objectiveStruct,r_exp_Sys,r_exp_Dias,S,VT_MRI,dispDias,dispSys)

%% Get variables
indVer=objectiveStruct.innerRingInd;
VT_img = objectiveStruct.VT_img;
H_opt = objectiveStruct.hOpt;
O_opt= objectiveStruct.oOpt;

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
V_unDEF = VT_MRI - K*U_find;

V_unDEF_innerRing = V_unDEF(indVer,:);
mean_V_unDEF = mean(V_unDEF_innerRing);
allRadialVecs_unDEF = V_unDEF_innerRing - (mean_V_unDEF).*ones (size(V_unDEF_innerRing));
meanRadius_unDEF = mean(sqrt(sum( allRadialVecs_unDEF.^2,2)));

%%  Inflate till diastole
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V_unDEF,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V_unDEF; %The nodel coordinates
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V_unDEF,1);

febio_spec.Material.material{1}.c1= C; %setting new stiffness
febio_spec.Material.material{1}.c2= C; %setting new stiffness

%Inflate till diastole
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureDias];

%Boundary conditions for diastole
febio_spec.Boundary.bc{2}.scale.VAL=dispDias(1);
febio_spec.Boundary.bc{3}.scale.VAL=dispDias(2);
febio_spec.Boundary.bc{4}.scale.VAL=dispDias(3);

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
febio_spec.Material.material{1}.c2= C; %setting new stiffness

febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;  1 appliedPressureSys];

%Boundary conditions for systole
febio_spec.Boundary.bc{2}.scale.VAL=dispSys(1);
febio_spec.Boundary.bc{3}.scale.VAL=dispSys(2);
febio_spec.Boundary.bc{4}.scale.VAL=dispSys(3);

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
    
        Popt = abs((r_exp_Sys - meanRadius_sys) - (r_exp_Dias - meanRadius_dias));
%     vDev= VT_img- V_DEF_dias;
%     Popt = (mean(sqrt(sum(vDev.^2,2))) + (r_exp_Sys - meanRadius_sys));
    
    OPT_stats_out.Popt=Popt;
    
    if ~isempty(H_opt)
        O_opt.XData= [O_opt.XData O_opt.XData(end)+1];
        O_opt.YData= [O_opt.YData (Popt)];
        H_opt.XData= [H_opt.XData H_opt.XData(end)+1];
        H_opt.YData= [H_opt.YData Pn];
        
        drawnow;
    end
    
    
    
else
    Popt=NaN(2,1); %Squared differences
    OPT_stats_out=[];
end
end



















































%% FUNCTIONS %%

%% Read lumen boundaries
function [outputStruct]=import_ctgr(fileName)

%% Access XML object
ctgrXML = xmlread(fileName);

%% Parse XML object

timestep_level = ctgrXML.getElementsByTagName('timestep');

numTimeSteps=timestep_level.getLength; %Number of time steps

for iTime=1:1:numTimeSteps
    contour_level = timestep_level.item(iTime-1).getElementsByTagName('contour');
    numContours=contour_level.getLength;
    for iContour=1:1:numContours
        contour_points_level = contour_level.item(iContour-1).getElementsByTagName('contour_points');
        point_level = contour_points_level.item(0).getElementsByTagName('point');
        num_point=point_level.getLength;
        
        V_points=zeros(num_point,3); %Allocate array for current point set
        for iPoint=1:1:num_point
            % id_points(iPoint)=str2double(point_level.item(iPoint-1).getAttribute('id').toCharArray()');
            V_points(iPoint,:)=[sscanf(point_level.item(iPoint-1).getAttribute('x').toCharArray()','%f')...
                sscanf(point_level.item(iPoint-1).getAttribute('y').toCharArray()','%f')...
                sscanf(point_level.item(iPoint-1).getAttribute('z').toCharArray()','%f')];
        end
        outputStruct.contour{iContour}.contour_points=V_points;
    end
end

end

%% Read center line points
function [outputStruct]=import_pth(fileName)

%% Access XML object
pthXML = xmlread(fileName);

%% Parse XML object

timestep_level = pthXML.getElementsByTagName('timestep');

numTimeSteps=timestep_level.getLength; %Number of time steps

for iTime=1:1:numTimeSteps
    controlpoints_level = timestep_level.item(iTime-1).getElementsByTagName('control_points');
    numControlpoints=controlpoints_level.getLength;
    for q=1:1:numControlpoints
        point_level = controlpoints_level.item(numControlpoints-1).getElementsByTagName('point');
        num_point=point_level.getLength;
        
        V_points=zeros(num_point,3); %Allocate array for current point set
        for iPoint=1:1:num_point
            % id_points(iPoint)=str2double(point_level.item(iPoint-1).getAttribute('id').toCharArray()');
            V_points(iPoint,:)=[sscanf(point_level.item(iPoint-1).getAttribute('x').toCharArray()','%f')...
                sscanf(point_level.item(iPoint-1).getAttribute('y').toCharArray()','%f')...
                sscanf(point_level.item(iPoint-1).getAttribute('z').toCharArray()','%f')];
        end
        outputStruct.centrePoints{q}.centroids=V_points;
    end
end

end


function [Vn]=curveOffset(V,wallThickness)
p1=mean(V,1); %Curve center

vf=-vecnormalize(V-[V(2:end,:);V(1,:)]); %Allong curve path vectors
vb=vecnormalize(V-[V(end,:);V(1:end-1,:)]); %Allong curve path vectors
v=(vf+vb)/2;


r=vecnormalize(V-p1); %Position vector wrt mean
v1=vecnormalize(cross(v,r)); %perimeter quasi Z-vectors
n=vecnormalize(cross(v1(ones(size(v,1),1),:),v)); %Outward normal vectors
Vn=(V+n*wallThickness); %Offset to create new curve

end

%% Resampling curve
function V=resampleCurve(V,pointSpacing,closeLoopOpt)

np=ceil(max(pathLength(V))/pointSpacing);
V = evenlySampleCurve(V,np,'spline',closeLoopOpt);

end




%%
function [Ve]=fitEllipse(V,n)

% Fit ellipse
[A] = ellipseFit3(V);

% Get coordinates for plotting
t=linspace(0,2*pi,n+1)'; t=t(1:end-1);
Ve=ellipseCoord3(A,t);

end

%%

function [F,V,CC,Cf]=cutMesh(Ft,Vt,Ct,Cf,Vc,Vc_ori,optionStructRayTrace,cutDist,bezierTangency,numSmootIterations,smoothIncludeSteps)

[~,~,Nt]=patchNormal(Ft,Vt);

pointSpacing=mean(patchEdgeLengths(Ft,Vt));
branchBaseOffset=cutDist;

%Ray trace curve to surface
p1=mean(Vc,1);
n1=vecnormalize(mean(Vc_ori,1)-p1);

V_cut_ellipse_ray_traced=Vc;
for q=1:1:size(V_cut_ellipse_ray_traced,1)
    pp=triSurfRayTrace(V_cut_ellipse_ray_traced(q,:),n1,[Ft(:,[1 2 3]); Ft(:,[3 4 1]);],Vt,optionStructRayTrace);
    if size(pp,1)>1
        [~,indMin]=minDist(V_cut_ellipse_ray_traced(q,:),pp);
        V_cut_ellipse_ray_traced(q,:)=pp(indMin,:);
    else
        V_cut_ellipse_ray_traced(q,:)=pp;
    end
end

%Cut hole
[~,indRemove]=minDist(V_cut_ellipse_ray_traced,Vt);
indRemove=unique(indRemove);
logicFacesSelect=~any(ismember(Ft,indRemove),2);

Ft_cut_away=Ft(~logicFacesSelect,:);
Cset=mean(Ct(~logicFacesSelect));

ns=mean(patchNormal(Ft_cut_away,Vt),1);

Ft_precut=Ft(logicFacesSelect,:);
Ct_precut=Ct(logicFacesSelect,:);
Cf_precut=Cf(logicFacesSelect,:);

optionStruct.outputType='label';
[G,~,groupSize]=tesgroup(Ft_precut,optionStruct);
[~,indLargestGroup]=max(groupSize);
Ft_cut=Ft_precut(G==indLargestGroup,:);
Ct_cut=Ct_precut(G==indLargestGroup);
Cf_cut=Cf_precut(G==indLargestGroup);

Et_boundary_cut=patchBoundary(Ft_cut);
Et_boundary=patchBoundary(Ft);
Et_cut_boundary=Et_boundary_cut(~any(ismember(Et_boundary_cut,Et_boundary),2),:);
indCutCurve=edgeListToCurve(Et_cut_boundary);

clear optionStruct
optionStruct.numSeeds=numel(indCutCurve); %Number of seeds
optionStruct.waitBarOn=0; %Turn on/off waitbar
Dt_cut=meshDistMarch(Ft_cut,Vt,indCutCurve,optionStruct);
logicVerticesFar=Dt_cut>cutDist;
logicKeep=all(logicVerticesFar(Ft_cut),2);
Ft_cut=Ft_cut(logicKeep,:);
Ct_cut=Ct_cut(logicKeep,:);
Cf_cut=Cf_cut(logicKeep,:);

%Loft transition
Et_boundary_cut=patchBoundary(Ft_cut);
Et_boundary=patchBoundary(Ft);
Et_cut_boundary=Et_boundary_cut(~any(ismember(Et_boundary_cut,Et_boundary),2),:);
indCutCurve=edgeListToCurve(Et_cut_boundary);
V_branch_curve_trunk=Vt(indCutCurve(1:end-1),:);


N_branch_curve_trunk=Nt(indCutCurve(1:end-1),:);

N_branch_curve_trunk=N_branch_curve_trunk +...
    [N_branch_curve_trunk(2:end,:); N_branch_curve_trunk(1,:)]+...
    [N_branch_curve_trunk(end,:); N_branch_curve_trunk(1:end-1,:);] ;
N_branch_curve_trunk=vecnormalize(N_branch_curve_trunk);


[~,indNearest]=minDist(V_branch_curve_trunk(1,:),V_cut_ellipse_ray_traced);
if indNearest>1
    V_cut_ellipse_ray_traced=[V_cut_ellipse_ray_traced(indNearest:end,:); V_cut_ellipse_ray_traced(1:indNearest-1,:)];
end
V_cut_ellipse_ray_traced=evenlySampleCurve(V_cut_ellipse_ray_traced,size(V_branch_curve_trunk,1),'pchip',1);

Vcm=Vc-p1+mean(V_cut_ellipse_ray_traced,1);

[~,indNearest]=minDist(V_branch_curve_trunk(1,:),Vcm);
if indNearest>1
    Vc=[Vc(indNearest:end,:); Vc(1:indNearest-1,:)];
end
Vc=evenlySampleCurve(Vc,size(V_branch_curve_trunk,1),'pchip',1);

branchAngle=180*acos(dot(-n1,ns))./pi;

V_branch_curve_1=V_cut_ellipse_ray_traced-n1.*(branchBaseOffset/cosd(90-branchAngle));

d1=sum(sqrt(sum((V_branch_curve_1-V_branch_curve_trunk).^2,2)));
d2=sum(sqrt(sum((flipud(V_branch_curve_1)-V_branch_curve_trunk).^2,2)));
if d2<d1
    V_branch_curve_1=flipud(V_branch_curve_1);
    Vc=flipud(Vc);
end

d=nan(size(V_branch_curve_1,1),1);
for q=1:1:size(V_branch_curve_1,1)
    if q>1
        V_test=[V_branch_curve_1(q:end,:);V_branch_curve_1(1:q-1,:)];
    else
        V_test=V_branch_curve_1;
    end
    d(q)=sum(sqrt(sum((V_test-V_branch_curve_trunk).^2,2)));
end
[~,indMin]=min(d);
if indMin>1
    V_branch_curve_1=[V_branch_curve_1(indMin:end,:);V_branch_curve_1(1:indMin-1,:)];
    Vc=[Vc(indMin:end,:);Vc(1:indMin-1,:)];
end

Vi=[V_branch_curve_trunk(end,:); V_branch_curve_trunk; V_branch_curve_trunk(1,:)];
di=pathLength(Vi);
di=di./max(di);

Vii=[V_branch_curve_1(end,:); V_branch_curve_1; V_branch_curve_1(1,:)];
dii=pathLength(Vii);
dii=dii./max(dii);

for q=1:size(Vi,2)
    V_branch_curve_1(:,q)=interp1(dii,Vii(:,q),di(2:end-1),'pchip');
end

% E1=vecnormalize(V_branch_curve_trunk- [V_branch_curve_trunk(2:end,:); V_branch_curve_trunk(1,:)]);
% E2=vecnormalize([V_branch_curve_trunk(end,:); V_branch_curve_trunk(1:end-1,:)]- V_branch_curve_trunk);
% E=vecnormalize(E1/2+E2/2);
% C1=vecnormalize(cross(N_branch_curve_trunk,E));
C1=vecnormalize(V_branch_curve_trunk-mean(V_branch_curve_trunk,1));
D1=vecnormalize(cross(C1,N_branch_curve_trunk));
M1=vecnormalize(cross(D1,N_branch_curve_trunk));
M2=-n1(ones(size(M1,1),1),:);

pointSpacingNow=mean(diff(pathLength(V_branch_curve_1)));%0.5*mean(diff(pathLength(V_branch_curve_1)))+0.5*mean(diff(pathLength(V_branch_curve_trunk)));
[Fb1,Vb1,X,Y,Z]=bezierLoft(V_branch_curve_trunk,V_branch_curve_1,M1,M2,pointSpacingNow,bezierTangency);


clear optionStruct
optionStruct.closeLoopOpt=1;
optionStruct.patchType='quad';
[Fb2,Vb2]=polyLoftLinear(V_branch_curve_1,Vc,optionStruct);
Fb2=fliplr(Fb2);

[F,V,C]=joinElementSets({Ft_cut,Fb1,Fb2},{Vt,Vb1,Vb2});
[F,V]=mergeVertices(F,V);
[F,V]=patchCleanUnused(F,V);
Cf=[Cf_cut; (max(Cf(:))+1)*ones(size(Fb1,1)+size(Fb2,1),1)];

CC=[Ct_cut; Cset*ones(size(Fb1,1)+size(Fb2,1),1)];

indSmoothRegion=F(C==2,:);
for q=1:1:smoothIncludeSteps
    indSmooth=F(any(ismember(F,indSmoothRegion),2),:);
    indSmoothRegion=indSmooth;
end
indRigid=F(~ismember(F,indSmooth));

clear optionStruct
optionStruct.Method='HC';
optionStruct.n=numSmootIterations;
optionStruct.RigidConstraints=indRigid;
V=patchSmooth(F,V,[],optionStruct);


end

%%

function [F,V,X,Y,Z]=bezierLoft(P1,P4,N1,N4,pointSpacing,f)

%Estimate required number of points
D12=sqrt(sum((P1-P4).^2,2));
numPoints=ceil(mean(D12)./pointSpacing);
if numPoints<2
    numPoints=2;
end

P2=P1+D12.*f.*N1;
P3=P4-D12.*f.*N4;

%Optimize required number of points by check bezier lengths
d=zeros(size(P1,1),1);
for q=1:1:size(P1,1)
    p=[P1(q,:); P2(q,:); P3(q,:); P4(q,:)]; %Control points
    V_bezier=bezierCurve(p,numPoints*2); %Compute bezier curve
    d(q)=max(pathLength(V_bezier));
end
numPoints=ceil(mean(d)./pointSpacing);
if numPoints<2
    numPoints=2;
end

%Compute final bezier curves
X=zeros(numPoints,size(P1,1));
Y=zeros(numPoints,size(P1,1));
Z=zeros(numPoints,size(P1,1));
for q=1:1:size(P1,1)
    p=[P1(q,:); P2(q,:); P3(q,:); P4(q,:)]; %Control points
    V_bezier=bezierCurve(p,numPoints*2); %Compute bezier curve
    V_bezier=evenlySampleCurve(V_bezier,numPoints,'pchip'); %resample evenly
    X(:,q)=V_bezier(:,1);
    Y(:,q)=V_bezier(:,2);
    Z(:,q)=V_bezier(:,3);
end

%Create quad patch data
[F,V] = surf2patch(X,Y,Z);
I=[(1:size(Z,1)-1)' (1:size(Z,1)-1)' (2:size(Z,1))' (2:size(Z,1))' ];
J=[size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1)];
F_sub=sub2ind(size(Z),I,J);
F=[F;F_sub];
F=fliplr(F);

end


