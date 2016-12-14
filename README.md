# MSGR - Matching Satellite and Ground Radar

Based on an IDL code from Rob Warren.

## Usage

You should only modify the config.ini file to match your own configuration and then run python matchvol.py in a terminal.

## Data

About the satellite data, here a copy of a mail from Rob Warren explaining how to get them:
The website where you can download TRMM and GPM data is https://storm.pps.eosdis.nasa.gov/storm/data/Service.jsp?serviceName=RestrictedOrder. You need to register before you can order data. One you've entered your pre-registered email address you can enter the details of your order. You want the following options:
*	Under 'Order type' select 'Standalone order'
*	Under 'Coincidence' select 'None or Satellite-Ground Validation Site'
*	Under 'Options' select 'Subset Geographic Area' then 'Subset Geographically' (leave 'Include only swaths with...' blank), and also select 'Parameter Subsetting'
*	Under 'Product Type' select '2AKu' under 'Algorithm' and check the box that comes up below. Note that this is for GPM. For TRMM you have to make two orders: one for 2A23 and one for 2A25, but if you're only working with GPM you don't need to worry about that.
*	Under 'Temporal Criteria' set the range of dates you want data for.
*	Under 'Special Area Of Interest' specify the limits of your domain (it should encompass the 150km range ring of your radar(s)). Give it a 'Location Alias'.
*	Under 'Parameter Subset' choose the following: 'dataQuality' from 'scanStatus', 'landSurfaceType' and 'flagPrecip' from 'PRE', 'flagBB', 'heightBB', 'widthBB', 'qualityBB', 'typePrecip', and 'qualityTypePrecip' from 'CSF', and 'zFactorCorrected' from 'SLV'. Provide an identifier, then select 'No' for 'Do you want to generate Read and Write routines for this subset', and set 'HDF' as the 'Output Data Format'.
*	Under 'Search Results' select all the files by clicking the top-most check box.
*	Under 'Script Type' select whichever you want: I used 'FTP URL' but you may prefer another.
Hitting submit should then get you what you want.
