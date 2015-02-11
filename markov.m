% Read the file
seq = fastaread('HMC21_NT_011515.txt');
seq = seq.Sequence;

% Remove nucleotides labelled as 'L'
seq = strrep(seq, 'L', '');

% Here, we initialize the tables and variables
p = .98;
q = .999;
chars = ['A' 'C' 'G' 'T'];
possible = ['A' 'C' 'G' 'T' 'a' 'c' 'g' 't'];
emissions = [1 0 0 0 
             0 1 0 0 
             0 0 1 0 
             0 0 0 1 
             1 0 0 0 
             0 1 0 0 
             0 0 1 0 
             0 0 0 1 
             ];

transMatrix = [(0.180*p)	(0.274*p)	(0.426*p)	(0.120*p)	(1-p)/4     (1-p)/4     (1-p)/4     (1-p)/4
              (0.171*p)	    (0.368*p)	(0.274*p)	(0.188*p)	(1-p)/4 	(1-p)/4     (1-p)/4     (1-p)/4
              (0.161*p)     (0.339*p)	(0.375*p)	(0.125*p)	(1-p)/4     (1-p)/4     (1-p)/4     (1-p)/4
              (0.079*p) 	(0.355*p)	(0.384*p)	(0.182*p)	(1-p)/4     (1-p)/4     (1-p)/4     (1-p)/4
              (1-q)/4       (1-q)/4     (1-q)/4     (1-q)/4 	(0.300*q)	(0.205*q)	(0.285*q)	(0.210*q)
              (1-q)/4       (1-q)/4     (1-q)/4     (1-q)/4     (0.322*q)	(0.298*q)	(0.078*q)	(0.302*q)
              (1-q)/4       (1-q)/4     (1-q)/4     (1-q)/4     (0.248*q)	(0.246*q)	(0.298*q)	(0.208*q)
              (1-q)/4       (1-q)/4     (1-q)/4     (1-q)/4     (0.177*q)	(0.239*q)	(0.292*q)	(0.292*q)];
                
tranMatLog = log(transMatrix);
emissionLog = log(emissions);

% Number of States and Characters
noOfStates = 8; 
m = 4;  

len = length(seq);

% Begin States and Characters
L = zeros(noOfStates,len);

for j = 1:4
    if seq(1) == chars(j)
        characterNum = j; 
        break;
    end 
end
    
for colIndex = 1:noOfStates
   L(colIndex,1) = log(1/noOfStates) + emissionLog(colIndex,characterNum);
end   		

% Find the scores for both possibilities
ptr = zeros(noOfStates,len); 
pi = zeros(1,len);
column = zeros(1,noOfStates);

for x = 2:len
    
    % Display the progress completed in percentage
    if mod(x,50000) == 1
        disp((x/len)*100);
    end
    
    for y = 1:noOfStates
        
        % Calculation of probability for each State
        for colIndex=1:noOfStates
            column(colIndex) = L(colIndex,x-1) + tranMatLog(colIndex,y); 
        end
    
        % Now, we can store the maximum probability and location 
        [L(y,x), ptr(y,x)] = max(column);
    
       % Next, we add the Character probability that is observed 
        for j = 1:4
            if seq(x) == chars(j)
                characterNum = j; 
                break;
            end 
        end
        L(y,x) = L(y,x) + emissionLog(y,characterNum);
    end 
end

% End the State
for colIndex = 1:noOfStates
   column(colIndex) = L(colIndex,len);
end
[maximumscore, pi(len)] = max(column);

% Find the most probable path for the sequence
for x = len:-1:2
   pi(x-1) = ptr(pi(x),x);
end

% Now, with the known matches gene regions, we can print out the Island locations and
% open the gene file location
fid = fopen('Chr21.txt','r');
tline = fgetl(fid);
i = 0;

% Read the file into an array
while (ischar(tline))
    i=i+1;
    [genes{i,1}, remain] = strtok(tline,'	');
    [genes{i,2}, remain] = strtok(remain,'	');
    [genes{i,3}, remain] = strtok(remain,'	');
    [genes{i,4}, remain] = strtok(remain,'	');   
    tline = fgetl(fid);
end

beginIteration = 0;
endIteration = 0;
geneName='';
begin = 43507093;
endReached = 1;
k = 0;

for i = 1:len    
    if pi(i)>0 & pi(i)<5 & beginIteration==0
        beginIteration = i;
    end
    if pi(i)>4 & beginIteration~=0
        endIteration = i;
    end
    
    if endIteration~=0 && beginIteration~=0
        if (endIteration - beginIteration)>=200
            beginIteration = begin+beginIteration;
            endIteration = begin+endIteration;
            
            geneLocation = str2num(cell2mat(genes(endReached,1)));
           
            if beginIteration > geneLocation
                if (geneLocation-beginIteration<500)
                    geneName = cell2mat(genes(endReached,4));
                else
                    geneName = '';
                end
                
                if endReached+1 < size(genes,1)
                    endReached = endReached + 1;
                else
                    geneName = '';
                end
            end
            
            k = k + 1;
            x = sprintf('CpG Island %i: \t %i bp \t (%i - %i) \t %s ', k, endIteration-beginIteration, beginIteration, endIteration, geneName);
            geneName = 0;
            disp(x);
        end
        
        beginIteration = 0;
        endIteration = 0;
    end    
end
     