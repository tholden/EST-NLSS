cd ..;
OldPath = path;
try
    RunEstimationText = fileread( 'RunEstimation.m' );
    RunEstimationText = regexprep( RunEstimationText, '^%.*$', '', 'lineanchors', 'dotexceptnewline' );
    RunEstimationFirstLine = regexp( RunEstimationText, 'function .*$', 'lineanchors', 'dotexceptnewline', 'once', 'match' );
    RunEstimationText = regexprep( RunEstimationText, 'function .*$', '', 'lineanchors', 'dotexceptnewline', 'once' );
    RunEstimationText = regexprep( RunEstimationText, '^[\r\n]*', '' );
    
    RunSmoothingText = fileread( 'RunSmoothing.m' );
    RunSmoothingText = regexprep( RunSmoothingText, '^%.*$', '', 'lineanchors', 'dotexceptnewline' );
    RunSmoothingFirstLine = regexp( RunSmoothingText, 'function .*$', 'lineanchors', 'dotexceptnewline', 'once', 'match' );
    RunSmoothingText = regexprep( RunSmoothingText, 'function .*$', '', 'lineanchors', 'dotexceptnewline', 'once' );
    RunSmoothingText = regexprep( RunSmoothingText, '^[\r\n]*', '' );

    ReadMeText = fileread( 'README.md' );
    ReadMeText = strrep( ReadMeText, '`', '' );
    ReadMeText = strrep( ReadMeText, '**', '' );
    ReadMeLines = strsplit( ReadMeText, { '\f', '\n', '\r', '\v' }, 'CollapseDelimiters', false );
    TxtFileID = fopen( 'ReadMe.txt', 'w' );
    RunEstimationFileID = fopen( 'RunEstimation.m', 'w' );
    RunSmoothingFileID = fopen( 'RunSmoothing.m', 'w' );
    fprintf( RunEstimationFileID, '%s\n', RunEstimationFirstLine );
    fprintf( RunSmoothingFileID, '%s\n', RunSmoothingFirstLine );
    
    for StrIdx = 1 : length( ReadMeLines )
        ReadMeLine = ReadMeLines{ StrIdx };
        SpaceString = regexp( ReadMeLine, '^\s*(\*\s*|\d+\.\s*)?', 'emptymatch', 'once', 'match' );
        SpaceLength = length( SpaceString );
        ReadMeLineWords = strsplit( ReadMeLine( SpaceLength+1:end ), { ' ', '\t' } );
        fprintf( TxtFileID, '%s', SpaceString );
        fprintf( RunEstimationFileID, '%% %s', SpaceString );
        fprintf( RunSmoothingFileID, '%% %s', SpaceString );
        LinePosition = SpaceLength;
        for WrdIdx = 1 : length( ReadMeLineWords )
            ReadMeLineWord = ReadMeLineWords{ WrdIdx };
            ReadMeLineWordLength = length( ReadMeLineWord );
            if LinePosition + 1 + ReadMeLineWordLength > 100 && SpaceLength + 1 + ReadMeLineWordLength <= 100
                SpaceString = regexprep( SpaceString, '\S', ' ' );
                fprintf( TxtFileID, '\n%s', SpaceString );
                fprintf( RunEstimationFileID, '\n%% %s', SpaceString );
                fprintf( RunSmoothingFileID, '\n%% %s', SpaceString );
                LinePosition = SpaceLength;
            end
            fprintf( TxtFileID, '%s ', ReadMeLineWord );
            fprintf( RunEstimationFileID, '%s ', ReadMeLineWord );
            fprintf( RunSmoothingFileID, '%s ', ReadMeLineWord );
            LinePosition = LinePosition + 1 + ReadMeLineWordLength;
        end
        fprintf( TxtFileID, '\n' );
        fprintf( RunEstimationFileID, '\n' );
        fprintf( RunSmoothingFileID, '\n' );
    end
    fprintf( RunEstimationFileID, '\n%s', RunEstimationText );
    fprintf( RunSmoothingFileID, '\n%s', RunSmoothingText );
    fclose( TxtFileID );
    fclose( RunEstimationFileID );
    fclose( RunSmoothingFileID );
        
    delete ReadMe.pdf;
    try
        addpath( 'C:\Program Files (x86)\Pandoc' );
    catch
    end
    !pandoc README.md -f markdown_github -o ReadMe.pdf -N --toc --wrap=none --latex-engine=xelatex -V papersize=A4 -V fontsize=10pt -V lang=en-GB -V documentclass=article -V margin-left=2.54cm -V margin-right=2.54cm -V margin-top=2.54cm -V margin-bottom=2.54cm -V mainfont=TeXGyrePagella -V sansfont=TeXGyreAdventor -V monofont=TeXGyreCursor -V links-as-notes -V colorlinks
catch Error
    disp( Error );
end
path( OldPath );
cd DevTools;
