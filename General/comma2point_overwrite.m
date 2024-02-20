function    comma2point_overwrite( filespec )
% replaces all occurences of comma (",") with point (".") in a text-file.
% Note that the file is overwritten, which is the price for high speed.
% see http://se.mathworks.com/matlabcentral/answers/57276-import-data-file-with-comma-decimal-point#answer_69283
% by Per Isakson
file    = memmapfile( filespec, 'writable', true );
comma   = uint8(',');
point   = uint8('.');
file.Data( transpose( file.Data==comma) ) = point;
end