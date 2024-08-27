function [gbest,gbestval, fitcount,RecordT,RESTART] = T_C9(fhd,D, Max_nfe, XRmin, XRmax, subfolderPath, subfolder, targetbest,varargin)
stm = RandStream('swb2712','Seed',sum(100*clock));
RandStream.setGlobalStream(stm);
RESTART = 0;
% ps_ini = 18*D;
ps_ini = round(25*log(D)*sqrt(D));
ps_min = 5;
ps = ps_ini;
profile on
if length(XRmin) == 1
    Rmin = repmat(XRmin,1,D);
    Rmax = repmat(XRmax,1,D);
end
VRmin = repmat(Rmin,ps,1);
VRmax = repmat(Rmax,ps,1);

fidvec = cell2mat(varargin);
fid = fidvec(1);
runid = fidvec(2);
WriteFlag =false;
if WriteFlag
    name = fullfile(subfolderPath,[subfolder,'_fid_',num2str(fid),'_',num2str(D),'D_',num2str(runid),'.dat']);
    fout = fopen(name,'wt');
end
tic;
%initialize;
pos = VRmin + (VRmax-VRmin).*rand(ps,D);
pastval = feval(fhd,pos',varargin{:});
fitcount = ps;
% copy_posval = posval;
[gbestval,gbestid] = min(pastval);
gbest = pos(gbestid,:);
if WriteFlag
    fprintf(fout,'%d\t%.15f\n',1,gbestval-targetbest(fid));
end

Cr = zeros(ps,1) + 0.9;
%     Cr0 = zeros(1,ps);
memory_size = 6; memory_order = 1;
memory_MUF = 0.5.*ones(memory_size,1); memory_MUCr = 0.9.*ones(memory_size,1);
Archfactor = 1.6;

% pbest_rate_max = 0.2; pbest_rate_min = 0.05;
pbest_rate = 0.15;

decayA = [];
T0 = 1.0;
% decayRate = T0/56;
A = [];
A_eval = [];
gen = 1;


Xmax = 100;Xmin = -100;
counter = zeros(ps,1);
Mean1 = mean(pos);
n = 50+D;
DI_ini = sqrt(sum(pdist2(Mean1,pos)));
t1 = repmat(0.9,ps,1);
% t2 = repmat(0.9,ps,1);
%     Fbest = repmat(0.5,ps,1);
%     Crbest = repmat(0.9,ps,1);
while fitcount < Max_nfe
    %       Distance = pdist2(gbest,pos);
    [~,indexSel] = sort(pastval);

    %% pbest_rate的数值也可以动态调整----无变化
    pNP = max(round(pbest_rate*ps),2);
    lenSel = max(ceil(pNP*rand(1,ps)),1);
    pbestIndex = zeros(ps,1);
    for idx = 1:ps
        pbestIndex(idx) = indexSel(lenSel(idx));
    end
    pbestB = pos(pbestIndex,:);

    if gen>1 && length(indexSel)>ps
        pos(indexSel(ps+1:length(indexSel))',:) = [];
        pastval(indexSel(ps+1:length(indexSel))') = [];
        MUF(indexSel(ps+1:length(indexSel))') = [];
        MUCr(indexSel(ps+1:length(indexSel))') = [];
        Cr(indexSel(ps+1:length(indexSel))') = [];
        sf(indexSel(ps+1:length(indexSel))') = [];
        counter(indexSel(ps+1:length(indexSel))') = [];
        VRmax = VRmax((1:ps)',:);
        VRmin = VRmin((1:ps)',:);
        if size(A,1)>round(Archfactor*ps)

            A = [A,decayA];
            A = A(1:round(Archfactor*ps),:);
            %             A_eval = A_eval(1:round(Archfactor*ps),:);
            decayA = A(:,end);
            A = A(:,1:(end-1));
        end
    end
    rndBase = randperm(ps)';
    psExt = ps + size(A,1);
    rndSeq1 = ceil(rand(ps,1)*psExt);
    for ii = 1:ps
        while rndBase(ii)==ii
            rndBase(ii)=ceil(rand()*ps);
        end
        while rndSeq1(ii)==rndBase(ii) || rndSeq1(ii)==ii
            rndSeq1(ii) = ceil(rand()*psExt);
        end
    end
    posx = [pos;A];
    posr = pos(rndBase,:);
    posxr = posx(rndSeq1,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       F/Cr生成公式       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    memory_rand_index1 = ceil(memory_size*rand(ps,1));
    %       memory_rand_index2 = ceil(memory_size*rand(ps,1));
    MUF = memory_MUF(memory_rand_index1);
    MUCr = memory_MUCr(memory_rand_index1);

    %%for generating crossover rate Cr
    rand1 = rand(ps,1);

    for ii = 1:ps
        if rand1(ii) < t1(ii)
            Cr(ii) = normrnd(MUCr(ii),0.1);
        end
        if find(MUCr(ii) == -1)
            Cr(ii) = 0;
        end
    end
    if fitcount < 0.4 * Max_nfe
        Cr = max(Cr,0.6);
    else
        Cr = min(Cr, 1);
        Cr = max(Cr, 0);
    end

    %% for generating scal factor F
    if (fitcount < Max_nfe * 0.22)
        sf = 2./sqrt(2).*pi^(-1/3).*(1-MUF.^2).*exp(-MUF.^2/2) + 0.15 * sin(pi * (rand(ps, 1) - 0.8));
        sf = min(sf, 0.6);
    else
        sf = MUF + 0.1 * tan(pi * (rand(ps, 1) - 0.5));
        %         for ii = 1:ps
        %             if rand2(ii) < t2(ii)
        %                 sf(ii,1) = randCauchy(MUF(ii),0.1);
        %             end
        %         end
        label_sf = find(sf <= 0);
        while ~ isempty(label_sf)
            sf(label_sf) = MUF(label_sf) + 0.1 * tan(pi * (rand(length(label_sf), 1) - 0.5));
            label_sf = find(sf <= 0);
        end


    end
    sf = min(sf,1);


    %% 能不能直接对应：换掉哪些维度会更好，而不是随机生成的点去比较自适应的Cr值。因为随机生成的不确定性太高会稀释Cr自适应选择带来的可解释性。
    label=zeros(ps,D);
    rndVal = rand(ps,D);
    onemat = zeros(ps,D);
    for ii = 1:ps
        label(ii,:) = rndVal(ii,:)<=Cr(ii);
        indexJ = ceil(rand()*D);
        onemat(ii,indexJ) = 1;
    end
    label = label|onemat;
    %% 初始的变异机制是从当前的种群开始，基于pbestB进行周围搜索-----能不能改成直接从pbestB开始 ，后面再次扰动
    offset = sf.*(pbestB-pos)+sf.*(posr-posxr);
    copy_pos = pos;
    pos = pos + offset;    %变异
    pos(~label) = copy_pos(~label);   %交叉
    pos = ((pos>=VRmin)&(pos<=VRmax)).*pos...
        +(pos<VRmin).*((VRmin+copy_pos).*rand(ps,D)/2) ...
        +(pos>VRmax).*((VRmax+copy_pos).*rand(ps,D)/2);

    dis = offset.*label;  % udiff

    posval = feval(fhd,pos',varargin{:});
    fitcount = fitcount + ps;


    %% 原有bin只包含对后代更好的标记，能不能加上对于劣等的评价，更充分的利用评估结果、
    %% 可以用语句：[fitness, I] = min([fitness, children_fitness], [], 2);
    %% I=1表示更新之后更好 ； I=2表示更新之后变差了
    [posval, I] = min([posval ; pastval], [], 1);
    pos(I == 2,:) = copy_pos(I == 2,:);   %选择
    copy_pos = [];
    %     bin = (pbestval > posval)'; %new better
    diff = (pastval-posval)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      至此交叉变异结束         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       时间戳       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 可以改一下这个A 变成对于bin=1的部分进行贮存
    A = [A;pos(I==1,:)];
    %     A_eval = [A_eval;posval(I==1)'];

    lengthAdd = size(pos(I==1,:),1);
    decayA = [decayA;zeros(lengthAdd,1)];
    decaystar = size(pos(I==1,:),1)/size(A,1);
    decayRate = tand(decaystar)/decaystar;
    decayA = decayA +decayRate;

    if(numel(decayA)>0)
        maxDecayA = max(decayA);
    else
        maxDecayA = 0;
    end

    if size(A,1)>round(Archfactor*ps) || maxDecayA >T0
        MergeA = [A,decayA];
        indexDecay = (decayA>T0);
        MergeA(indexDecay,:) = [];
        %         A_eval(indexDecay,:) = [];
        len = length(MergeA(:,1));
        if len>round(Archfactor*ps)
            rndSel = randperm(len)';
            rndSel = rndSel(round(Archfactor*ps)+1:len);
            MergeA(rndSel,:) = [];
        end
        %         if len>round(Archfactor*ps)
        %             [A_eval,A_label] = sort(A_eval,"descend");
        %
        %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    MergeA = MergeA(A_label(1:round(Archfactor*ps)),:);
        %             A_eval = A_eval(A_label(1:round(Archfactor*ps)),1);
        %         end
        A = MergeA(:,1:D);
        decayA = MergeA(:,D+1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       时间戳结束       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       F\Cr参数自适应更新      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     pbest(bin==1,:) = pos(bin==1,:);
    %     pbestval(bin==1) = posval(bin==1);

    %         Fbest(bin==1) = sf(bin==1);
    %         Crbest(bin==1) = Cr(bin==1);

    SuccF = sf(I==1);
    SuccCr = Cr(I==1);
    %     diffB = diff(I==1);
    %             diffB = diffB';
    %     diffB = diffB/sum(diffB);
    dis = dis(I==1,:);   % usucc
    %     dis = std(dis')';
    %     dis = dis/sum(dis);
    FCR = std(dis,0,2);
    FCR = FCR/sum(FCR);

    %         Improve_num = label(I==1,:);
    %         Improve_num1 = sum(Improve_num,2)/D;
    %         label2 = (Improve_num1/sum(Improve_num1));
    %         label2 = label2';
    num_Succ = numel(SuccCr);
    if num_Succ > 0
        memory_MUF(memory_order) = (FCR'*(SuccF.^2))/(FCR'*SuccF);
        %              if fitcount < pivot(1)
        %memory_MUF(memory_order) = (dis'*(SuccF.^2))/(dis'*SuccF);
        %memory_MUF(memory_order) = (diffB'*(SuccF.^2))/(diffB'*SuccF);
        %memory_MUCr(memory_order) = (label2*(SuccF.^2))./(label2*SuccF);
        if max(SuccCr) == 0 || memory_MUCr(memory_order) == -1
            memory_MUCr(memory_order) = -1;
        else
            memory_MUCr(memory_order) = ((SuccCr.^2)'*FCR)/(SuccCr'*FCR);
            %             memory_MUCr(memory_order) = (diffB'*(SuccCr.^2))/(diffB'*SuccCr);
            %             memory_MUCr(memory_order) = (label2*(SuccCr.^2))./(label2*SuccCr);
            %memory_MUCr(memory_order) = (dis'*(SuccCr.^2))/(dis'*SuccCr);
        end

        memory_order = memory_order + 1;
        if memory_order > memory_size
            memory_order = 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       F\Cr参数自适应更新结束      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [gbestval,gbestid] = min(posval);
    gbest = pos(gbestid,:);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      重启机制      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %% restart mechanism
     for ii = 1:ps
       if I(ii) == 2
            counter(ii) = counter(ii) + 1;
        else
            counter(ii) = 0;
        end
    end
    V_STD = std(pos);
    [~,v_label] = sort(V_STD,'descend');
    v_label = v_label(1:ceil(0.3*D));
    Meanpos = mean(pos);
    DI = sqrt(sum(pdist2(Meanpos,pos)));
    RDI = DI/DI_ini;
    if RDI < 0.02
        for i = 1:ps
            if counter(i)>=n && i~=gbestid(1)
                pos(i,v_label) = gbest(v_label) + 0.1 *rand(1,length(v_label)).* ( gbest(v_label) - Meanpos(v_label) );
                posval(i) = feval(fhd,pos(i,:)',fid);
                counter(i)=0;
                fitcount=fitcount+1;
                RESTART = RESTART+1;
            end
        end
    end
%        C_o = 0;
%     for ii = 1:ps
%         if I(ii) == 2
%             counter(ii) = counter(ii) + 1;
%             C_o = C_o + counter(ii);
%         else
%             counter(ii) = 0;
%         end
%     end
%   V_STD = std(pos);
%     [~,v_label] = sort(V_STD,'descend');
%     v_label = v_label(1:8);
%     Meanpos = mean(pos);
%     DI = sqrt(sum(pdist2(Meanpos,pos)));
%     RDI = DI/DI_ini;
%     if RDI < 0.001
%         for i = 1:ps
%             if counter(i)>=n && i~=gbestid(1)  &&  C_o >0.7 * ps * D%C_o > num_Co其中num_Co = 0.7 * ps * D;
%                 %             rand_num = ceil(ps/10);
%                 pos (i,v_label) = gbest(v_label) + 0.1 *rand(1,length(v_label)).* ( gbest(v_label) - Meanpos(v_label) );
%                 posval(i) = feval(fhd,pos(i,:)',fid);
%                 counter(i)=0;
%                 fitcount=fitcount+1;
%                 RESTART = RESTART+1;
%             end
%         end
%         C_o = 0;
%     end



%     C_o = 0;
%     for i = 1:ps
%         if I(i) == 2
%             counter(i) = counter(i) + 1;
%             C_o = C_o + counter(i);
%         else
%             counter(i) = 0;
%         end
%     end
%     for v=1:D
%         V_upon(v) = max(pos(:,v));
%         V_down(v) = min(pos(:,v));
%     end
%     V_pos = V_upon -V_down;
%     V_num = 0;
%     for j=1:D
%         if V_pos(j) < 0.1
%             V_num =V_num + 1;
%         end
%     end
%     if V_num == ps
%         for i = 1:ps
%             if  C_o > num_Co  && i~=gbestid && counter(i) > n
% 
%                 nDim = randi(ps);
%                 nSeq=randperm(ps);
%                 j=nSeq(1:nDim);
%                 [~,j_size] = size(j);
%                 for j1 = 1:j_size
%                     j_num = j(1,j1);
%                     pop_num = pos(:,j_num);
%                     pop_num_1 = pop_num(randperm(numel(pop_num),1));
%                     popold(i,j_num)=pop_num_1;
% 
%                 end
%                 fitness(i) = feval(fhd,popold(i,:)',func);
%                 counter(i)=0;
%                 fitcount = fitcount + 1;
%             end
% 
%         end
%     end

    %     Meanpos = mean(pos);
    %     DI = sqrt(sum(pdist2(Meanpos,pos)));
    %     RDI = DI/DI_ini;
    %     if RDI < 0.001
    %         for i = 1:ps
    %             if counter(i)>=n && i~=gbestid
    %                 %% 可以改： 重启策略可以基于种群聚集区域生成，这样可以更好的服务于勘探
    %                 %% 原来是测试向量重新生成，现在改成维度重排序策略--感觉没什么道理，先试试
    %                 %%  C_o = C_o + counter(i);  C_o >  n  && i~=gbestid && counter(i) > 2*problem_size
    %                 nDim = randi(D);
    %                 sequence = randperm(D, nDim); %随便抽取了一段维度序列进行排序
    %                 copy_pos = pos;
    %                 for j= 1: nDim
    %                     CrossingObject = randi(ps);
    %                     pos(i,j) = copy_pos(CrossingObject,j);
    %                 end
    %                 clear copy_pos;
    %
    %                 %% 注意这里设计大量评估次数，是不是可以优化一下
    %                 posval(i) = feval(fhd,pos(i,:)',varargin{:});%pastval = feval(fhd,pos',varargin{:});
    %                 counter(i)=0;
    %                 fitcount=fitcount+1;
    %             end
    %         end
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      重启机制结束      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        if((mod(fitcount-ps,100*5*D)>=100*5*D-ps)&& mod(fitcount,100*5*D)<ps)  &&  WriteFlag
            %if mod(gen,5*D)==0
            fprintf(fout,'%d\t%.15f\n',fitcount,gbestval-targetbest(fid));
        end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      种群递减策略      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         pivot = [1/3*Max_nfe,2/3*ps_ini];
% 
%         if fitcount< pivot(1)
%             plan_ps = ceil((pivot(2)-ps_ini)/(pivot(1)-ps_ini)^2*(fitcount-ps_ini)^2 + ps_ini);
%         else
%             plan_ps = floor((pivot(2)-ps_min)/(pivot(1)-Max_nfe)*(fitcount-Max_nfe)+ps_min);
%         end
    pivot = [2/3*Max_nfe,1/3*ps_ini];

    if fitcount< ceil(pivot(1))
        plan_ps = ceil((pivot(2)-ps_ini)/(pivot(1)-ps_ini)^2*(fitcount-ps_ini)^2 + ps_ini);
    else
        plan_ps = floor((pivot(2)-ps_ini)/(pivot(1)-ps_ini)*(fitcount-Max_nfe)+ps_min);
    end
    % plan_ps = round((((ps_min - ps_ini) / Max_nfe) * fitcount) + ps_ini);
    reduc_nums = ps - plan_ps;
    %% 记录ps变化情况
    %         reduce_ps = fullfile(subfolderPath,['TEST03_', num2str(fff), '_reduce-ps.txt']);
    %         f_reduce = fopen(reduce_ps,'a');
    %        fprintf(f_reduce,'%d\n',reduc_nums);

    if ps > plan_ps
        if ps - reduc_nums < ps_min
            reduc_nums = ps - ps_min;
        end
        ps = ps - reduc_nums;
    end
    %     pbest_rate = pbest_rate_max+(pbest_rate_min-pbest_rate_max)*fitcount/Max_nfe;
    pbest_rate =  0.15 + 0.3*fitcount/Max_nfe;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      种群递减策略结束      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gen = gen +1;
    pastval = posval;
end
RecordT = toc;

if WriteFlag
    fclose(fout);
end
profile off%profile viewer
end

function result = randCauchy(mu, sigma)
[m,n] = size(mu);
result = mu + sigma*tan(pi*(rand(m,n)-0.5));
end