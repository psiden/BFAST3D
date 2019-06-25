classdef HutchRandomPool < handle
  %HutchRandomPool A pool of random vectors Vs to be used for Hutchinson
  % estimation. Stores old solutions, for good starting guesses
  %
  % AUTHOR:       Per Siden
  %               Division of Statistics and Machine Learning
  %               Department of Computer and Information Science
  %               Linkoping University
  %
  % FIRST VER.:
  % REVISED.:     2017-09-15
  
  properties
    
    N                   % Vector length
    S                   % Number of samples in the pool
    Ns                  % Number of samples to be used in Hutch
    Vs                  % Random vectors
    lastiQVs            % Last solutions
    lastInd             % Sample indices used
    
  end
  
  methods
    
    function ob = HutchRandomPool(N,S,Ns)
      assert(S>=Ns);
      ob.N = N;
      ob.S = S;
      ob.Ns = Ns;
      ob.Vs = 2*(round(rand(N,S))-.5);
      ob.lastiQVs = zeros(N,S);
      
    end
    
    function [Vs,lastiQVs] = getSamples(ob)
      ob.lastInd = randsample(ob.S,ob.Ns,false);
      Vs = ob.Vs(:,ob.lastInd);
      lastiQVs = ob.lastiQVs(:,ob.lastInd);
    end
    
    function storeSolutions(ob,iQVs)
      ob.lastiQVs(:,ob.lastInd) = iQVs;
    end
    
  end
end

