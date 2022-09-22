function se=strel3d(sesize)

% function se=STREL3D(sesize)
%
% STREL3D creates a 3D sphere as a structuring element. Three-dimensional 
% structuring elements are much better for morphological reconstruction and
% operations of 3D datasets. Otherwise the traditional MATLAB "strel"
% function will only operate on a slice-by-slice approach. This function
% uses the aribtrary neighborhood for "strel."
% 
% Usage:        se=STREL3D(sesize)
%
% Arguments:    sesize - desired diameter size of a sphere (any positive 
%               integer)
%
% Returns:      the structuring element as a strel class (can be used
%               directly for imopen, imclose, imerode, etc)
% 
% Examples:     se=strel3d(1)
%               se=strel3d(2)
%               se=strel3d(5)
%
% 2014/09/26 - LX 
% 2014/09/27 - simplification by Jan Simon


sw=(sesize-1)/2; 
ses2=ceil(sesize/2);            % ceil sesize to handle odd diameters
[y,x,z]=meshgrid(-sw:sw,-sw:sw,-sw:sw); 
m=sqrt(x.^2 + y.^2 + z.^2); 
b=(m <= m(ses2,ses2,sesize)); 
se=strel('arbitrary',b);



% citation
% Luke Xie (2022). 3D structuring element (sphere) (https://www.mathworks.com/matlabcentral/fileexchange/47937-3d-structuring-element-sphere), MATLAB Central File Exchange

% License
% Copyright (c) 2014, Luke Xie
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of  nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.