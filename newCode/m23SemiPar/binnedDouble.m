classdef binnedDouble < double
    properties (SetAccess = public)
        IX={};
    end
   methods
      function obj = binnedDouble(x,IX)
         obj = obj@double(x);
         IX=IX(:);
         for i=1:length(IX)
             IX{i}=IX{i}(:);
         end
         obj.IX = IX;  
      end

       function obj = concat(obj1,obj2)
          n=length(obj1);
          x= [double(obj1);double(obj2)];
          IX=cell(length(obj1.IX),1);
          for i= 1:length(obj1.IX)
              IX{i}=[obj1.IX{i};obj2.IX{i}+n];
          end
          obj = binnedDouble(x,IX);
       end
       function [x,ix] = sort(obj)
          [x,ix]= sort(double(obj));
          n=0;
          IX=cell(length(obj.IX),1);
          for i= 1:length(obj.IX)
              ni=length(obj.IX{i});
              IX{i}=n+ (1:ni);
              n = n + ni;
          end
          x = binnedDouble(x,IX);
      end
      
   end
end