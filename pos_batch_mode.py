import os


def generate_positive_batch_mode(
    file_name: str,
    feature_list: str,
    output_sirius: str,
    output_gnps: str,
) -> str:

    return f"""
<?xml version="1.0" encoding="UTF-8"?><batch mzmine_version="4.6.1">
    <batchstep method="io.github.mzmine.modules.io.import_rawdata_all.AllSpectralDataImportModule" parameter_version="1">
        <parameter name="File names">
            <file>{file_name}</file>
        </parameter>
        <parameter name="Try vendor centroiding">true</parameter>
        <parameter name="Advanced import" selected="false">
            <parameter name="Scan filters" selected="true">
                <parameter name="Scan number"/>
                <parameter name="Base Filtering Integer"/>
                <parameter name="Retention time"/>
                <parameter name="Mobility"/>
                <parameter name="MS level filter" selected="All MS levels">1</parameter>
                <parameter name="Scan definition"/>
                <parameter name="Polarity">Any</parameter>
                <parameter name="Spectrum type">ANY</parameter>
            </parameter>
            <parameter name="Crop MS1 m/z" selected="false"/>
            <parameter name="MS1 detector (Advanced)" selected="false" selected_item="Factor of lowest signal">
                <module name="Factor of lowest signal">
                    <parameter name="Noise factor">2.5</parameter>
                </module>
                <module name="Auto">
                    <parameter name="Noise level">1000.0</parameter>
                </module>
                <module name="Centroid">
                    <parameter name="Noise level"/>
                </module>
                <module name="Exact mass">
                    <parameter name="Noise level"/>
                </module>
                <module name="Local maxima">
                    <parameter name="Noise level"/>
                </module>
                <module name="Recursive threshold">
                    <parameter name="Noise level"/>
                    <parameter name="Min m/z peak width"/>
                    <parameter name="Max m/z peak width"/>
                </module>
                <module name="Wavelet transform">
                    <parameter name="Noise level"/>
                    <parameter name="Scale level"/>
                    <parameter name="Wavelet window size (%)"/>
                </module>
            </parameter>
            <parameter name="MS2 detector (Advanced)" selected="false" selected_item="Factor of lowest signal">
                <module name="Factor of lowest signal">
                    <parameter name="Noise factor">2.5</parameter>
                </module>
                <module name="Auto">
                    <parameter name="Noise level">1000.0</parameter>
                </module>
                <module name="Centroid">
                    <parameter name="Noise level"/>
                </module>
                <module name="Exact mass">
                    <parameter name="Noise level"/>
                </module>
                <module name="Local maxima">
                    <parameter name="Noise level"/>
                </module>
                <module name="Recursive threshold">
                    <parameter name="Noise level"/>
                    <parameter name="Min m/z peak width"/>
                    <parameter name="Max m/z peak width"/>
                </module>
                <module name="Wavelet transform">
                    <parameter name="Noise level"/>
                    <parameter name="Scale level"/>
                    <parameter name="Wavelet window size (%)"/>
                </module>
            </parameter>
            <parameter name="Denormalize fragment scans (traps)">false</parameter>
        </parameter>
        <parameter name="Metadata file" selected="false"/>
        <parameter name="Sort and color">true</parameter>
        <parameter name="Spectral library files"/>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetectionModule" parameter_version="1">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scan filters" selected="true">
            <parameter name="Scan number"/>
            <parameter name="Base Filtering Integer"/>
            <parameter name="Retention time">
                <min>0.0</min>
                <max>12.01</max>
            </parameter>
            <parameter name="Mobility"/>
            <parameter name="MS level filter" selected="MS2, level = 2">1</parameter>
            <parameter name="Scan definition"/>
            <parameter name="Polarity">+</parameter>
            <parameter name="Spectrum type">ANY</parameter>
        </parameter>
        <parameter name="Scan types (IMS)">All scan types</parameter>
        <parameter name="Denormalize fragment scans (traps)">false</parameter>
        <parameter name="Mass detector" selected_item="Centroid">
            <module name="Factor of lowest signal">
                <parameter name="Noise factor">2.5</parameter>
            </module>
            <module name="Auto">
                <parameter name="Noise level">1000.0</parameter>
            </module>
            <module name="Centroid">
                <parameter name="Noise level">0.0</parameter>
            </module>
            <module name="Exact mass">
                <parameter name="Noise level"/>
            </module>
            <module name="Local maxima">
                <parameter name="Noise level"/>
            </module>
            <module name="Recursive threshold">
                <parameter name="Noise level"/>
                <parameter name="Min m/z peak width"/>
                <parameter name="Max m/z peak width"/>
            </module>
            <module name="Wavelet transform">
                <parameter name="Noise level"/>
                <parameter name="Scale level"/>
                <parameter name="Wavelet window size (%)"/>
            </module>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_massdetection.MassDetectionModule" parameter_version="1">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scan filters" selected="true">
            <parameter name="Scan number"/>
            <parameter name="Base Filtering Integer"/>
            <parameter name="Retention time">
                <min>0.0</min>
                <max>12.01</max>
            </parameter>
            <parameter name="Mobility"/>
            <parameter name="MS level filter" selected="MS1, level = 1">1</parameter>
            <parameter name="Scan definition"/>
            <parameter name="Polarity">+</parameter>
            <parameter name="Spectrum type">ANY</parameter>
        </parameter>
        <parameter name="Scan types (IMS)">All scan types</parameter>
        <parameter name="Denormalize fragment scans (traps)">false</parameter>
        <parameter name="Mass detector" selected_item="Centroid">
            <module name="Factor of lowest signal">
                <parameter name="Noise factor">2.5</parameter>
            </module>
            <module name="Auto">
                <parameter name="Noise level">1000.0</parameter>
            </module>
            <module name="Centroid">
                <parameter name="Noise level">10000.0</parameter>
            </module>
            <module name="Exact mass">
                <parameter name="Noise level"/>
            </module>
            <module name="Local maxima">
                <parameter name="Noise level"/>
            </module>
            <module name="Recursive threshold">
                <parameter name="Noise level"/>
                <parameter name="Min m/z peak width"/>
                <parameter name="Max m/z peak width"/>
            </module>
            <module name="Wavelet transform">
                <parameter name="Noise level"/>
                <parameter name="Scale level"/>
                <parameter name="Wavelet window size (%)"/>
            </module>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.featdet_targeted.TargetedFeatureDetectionModule" parameter_version="1">
        <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
        <parameter name="Scan filters" selected="true">
            <parameter name="Scan number"/>
            <parameter name="Base Filtering Integer"/>
            <parameter name="Retention time">
                <min>0.0</min>
                <max>12.01</max>
            </parameter>
            <parameter name="Mobility"/>
            <parameter name="MS level filter" selected="MS1, level = 1">1</parameter>
            <parameter name="Scan definition"/>
            <parameter name="Polarity">+</parameter>
            <parameter name="Spectrum type">ANY</parameter>
        </parameter>
        <parameter name="Name suffix">detectedPeak</parameter>
        <parameter name="Database file">
            <current_file>{feature_list}</current_file>
        </parameter>
        <parameter name="Field separator">,</parameter>
        <parameter name="Columns">
            <importtype column="neutral mass" datatype="io.github.mzmine.datamodel.features.types.numbers.NeutralMassType" selected="true"/>
            <importtype column="mz" datatype="io.github.mzmine.datamodel.features.types.numbers.PrecursorMZType" selected="false"/>
            <importtype column="rt" datatype="io.github.mzmine.datamodel.features.types.numbers.RTType" selected="true"/>
            <importtype column="formula" datatype="io.github.mzmine.datamodel.features.types.annotations.formula.FormulaType" selected="true"/>
            <importtype column="smiles" datatype="io.github.mzmine.datamodel.features.types.annotations.SmilesStructureType" selected="false"/>
            <importtype column="adduct" datatype="io.github.mzmine.datamodel.features.types.annotations.iin.IonAdductType" selected="false"/>
            <importtype column="inchi" datatype="io.github.mzmine.datamodel.features.types.annotations.InChIStructureType" selected="false"/>
            <importtype column="inchi key" datatype="io.github.mzmine.datamodel.features.types.annotations.InChIKeyStructureType" selected="false"/>
            <importtype column="name" datatype="io.github.mzmine.datamodel.features.types.annotations.CompoundNameType" selected="true"/>
            <importtype column="CCS" datatype="io.github.mzmine.datamodel.features.types.numbers.CCSType" selected="false"/>
            <importtype column="mobility" datatype="io.github.mzmine.datamodel.features.types.numbers.MobilityType" selected="false"/>
            <importtype column="comment" datatype="io.github.mzmine.datamodel.features.types.annotations.CommentType" selected="false"/>
        </parameter>
        <parameter name="Intensity tolerance">0.1</parameter>
        <parameter name="m/z tolerance (scan-to-scan)">
            <absolutetolerance>0.01</absolutetolerance>
            <ppmtolerance>10.0</ppmtolerance>
        </parameter>
        <parameter name="Retention time tolerance" selected="true" unit="MINUTES">12.0</parameter>
        <parameter name="Mobility time tolerance" selected="false"/>
        <parameter name="Calculate adduct masses" selected="true">
            <parameter name="Maximum charge">2</parameter>
            <parameter name="Maximum molecules/cluster">3</parameter>
            <parameter name="Adducts">
                <adduct_type selected="false">
                    <subpart charge="1" mass_difference="-5.4858E-4" mol_formula="" name="e" type="ADDUCT"/>
                    <subpart charge="1" mass_difference="-5.4858E-4" mol_formula="" name="e" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="-1" mass_difference="-1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="1" mass_difference="22.989218" mol_formula="Na" name="Na" type="ADDUCT"/>
                    <subpart charge="1" mass_difference="22.989218" mol_formula="Na" name="Na" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="0" mass_difference="-18.010565" mol_formula="H2O" name="H2O" type="NEUTRAL_LOSS"/>
                    <subpart charge="1" mass_difference="-5.4858E-4" mol_formula="" name="e" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="0" mass_difference="-18.010565" mol_formula="H2O" name="H2O" type="NEUTRAL_LOSS"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="2" mass_difference="55.93384" mol_formula="Fe" name="Fe" type="ADDUCT"/>
                    <subpart charge="-1" mass_difference="-1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="1" mass_difference="18.033823" mol_formula="NH4" name="NH4" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="2" mass_difference="39.96149382" mol_formula="Ca" name="Ca" type="ADDUCT"/>
                    <subpart charge="-1" mass_difference="-1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="1" mass_difference="38.963158" mol_formula="K" name="K" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="2" mass_difference="47.96953482" mol_formula="Mg" name="Mg" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="1" mass_difference="38.963158" mol_formula="K" name="K" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="1" mass_difference="22.989218" mol_formula="Na" name="Na" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="1" mass_difference="-5.4858E-4" mol_formula="" name="e" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="1" mass_difference="22.989218" mol_formula="Na" name="Na" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="2" mass_difference="55.93384" mol_formula="Fe" name="Fe" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="0" mass_difference="-18.010565" mol_formula="H2O" name="H2O" type="NEUTRAL_LOSS"/>
                    <subpart charge="0" mass_difference="-18.010565" mol_formula="H2O" name="H2O" type="NEUTRAL_LOSS"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="2" mass_difference="39.96149382" mol_formula="Ca" name="Ca" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="0" mass_difference="-18.010565" mol_formula="H2O" name="H2O" type="NEUTRAL_LOSS"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="0" mass_difference="-18.010565" mol_formula="H2O" name="H2O" type="NEUTRAL_LOSS"/>
                    <subpart charge="0" mass_difference="-18.010565" mol_formula="H2O" name="H2O" type="NEUTRAL_LOSS"/>
                    <subpart charge="0" mass_difference="-18.010565" mol_formula="H2O" name="H2O" type="NEUTRAL_LOSS"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="-1" mass_difference="-1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                    <subpart charge="2" mass_difference="47.96953482" mol_formula="Mg" name="Mg" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="1" mass_difference="18.033823" mol_formula="NH4" name="NH4" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="true">
                    <subpart charge="1" mass_difference="1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="-1" mass_difference="-1.007276" mol_formula="H" name="H" type="ADDUCT"/>
                </adduct_type>
                <adduct_type selected="false">
                    <subpart charge="1" mass_difference="0.0" mol_formula="" name="e" type="ADDUCT"/>
                </adduct_type>
            </parameter>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.dataprocessing.filter_groupms2.GroupMS2Module" parameter_version="3">
        <parameter name="Feature lists" type="BATCH_LAST_FEATURELISTS"/>
        <parameter name="MS1 to MS2 precursor tolerance (m/z)">
            <absolutetolerance>0.01</absolutetolerance>
            <ppmtolerance>10.0</ppmtolerance>
        </parameter>
        <parameter name="Retention time filter" selected="Use feature edges" unit="MINUTES">0.2</parameter>
        <parameter name="Minimum relative feature height" selected="true">0.25</parameter>
        <parameter name="Minimum required signals" selected="true">1</parameter>
        <parameter name="Limit by ion mobility edges">false</parameter>
        <parameter name="Minimum detections in IMS dimension">2</parameter>
        <parameter name="Merge MS/MS spectra (TIMS)">false</parameter>
        <parameter name="Advanced" selected="false">
            <parameter name="Minimum signal intensity (absolute, TIMS)" selected="false">250.0</parameter>
            <parameter name="Minimum signal intensity (relative, TIMS)" selected="true">0.01</parameter>
            <parameter name="Group iterative MS2s" selected="false">mainQuantFile</parameter>
        </parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.io.export_features_sirius.SiriusExportModule" parameter_version="2">
        <parameter name="Feature lists" type="BATCH_LAST_FEATURELISTS"/>
        <parameter name="Filename">
            <current_file>{output_sirius}</current_file>
        </parameter>
        <parameter name="Intensity normalization" scientific="true">no_normalization</parameter>
        <parameter name="Merge &amp; select fragment scans" selected_item="simple_merged">
            <module name="simple_merged">
                <parameter name="Presets">single_merged_scan</parameter>
                <parameter name="Merging m/z tolerance">
                    <absolutetolerance>0.001</absolutetolerance>
                    <ppmtolerance>5.0</ppmtolerance>
                </parameter>
            </module>
            <module name="preset_merged">
                <parameter name="Presets">representative_scans</parameter>
                <parameter name="Merging m/z tolerance">
                    <absolutetolerance>0.008</absolutetolerance>
                    <ppmtolerance>25.0</ppmtolerance>
                </parameter>
                <parameter name="Merge">
                    <selected>Across samples</selected>
                </parameter>
                <parameter name="Intensity merge mode">MAXIMUM</parameter>
            </module>
            <module name="input_scans">
                <parameter name="Select input scans">most_intense_across_samples</parameter>
            </module>
            <module name="advanced">
                <parameter name="Merging options">
                    <selected>Across samples</selected>
                    <selected>Across energies</selected>
                </parameter>
                <parameter name="m/z tolerance">
                    <absolutetolerance>0.008</absolutetolerance>
                    <ppmtolerance>25.0</ppmtolerance>
                </parameter>
                <parameter name="Intensity merge mode">MAXIMUM</parameter>
                <parameter name="Also include input scans">none</parameter>
            </module>
        </parameter>
        <parameter name="m/z tolerance">
            <absolutetolerance>0.003</absolutetolerance>
            <ppmtolerance>5.0</ppmtolerance>
        </parameter>
        <parameter name="Only rows with annotation">false</parameter>
        <parameter name="Exclude multiple charge">false</parameter>
        <parameter name="Exclude multimers">false</parameter>
    </batchstep>
    <batchstep method="io.github.mzmine.modules.io.export_features_gnps.fbmn.GnpsFbmnExportAndSubmitModule" parameter_version="3">
        <parameter name="Feature lists" type="BATCH_LAST_FEATURELISTS"/>
        <parameter name="Filename">
            <current_file>{output_gnps}</current_file>
        </parameter>
        <parameter name="Filter rows">ONLY WITH MS2</parameter>
        <parameter name="Merge &amp; select fragment scans" selected_item="input_scans">
            <module name="simple_merged">
                <parameter name="Presets">single_merged_scan</parameter>
                <parameter name="Merging m/z tolerance">
                    <absolutetolerance>0.008</absolutetolerance>
                    <ppmtolerance>25.0</ppmtolerance>
                </parameter>
            </module>
            <module name="preset_merged">
                <parameter name="Presets">single_merged_scan</parameter>
                <parameter name="Merging m/z tolerance">
                    <absolutetolerance>0.008</absolutetolerance>
                    <ppmtolerance>25.0</ppmtolerance>
                </parameter>
                <parameter name="Merge">
                    <selected>Across samples</selected>
                </parameter>
                <parameter name="Intensity merge mode">MAXIMUM</parameter>
            </module>
            <module name="input_scans">
                <parameter name="Select input scans">most_intense_across_samples</parameter>
            </module>
        </parameter>
        <parameter name="Intensity normalization" scientific="false">no_normalization</parameter>
        <parameter name="Feature intensity">Height</parameter>
        <parameter name="CSV export">SIMPLE</parameter>
        <parameter name="Submit to GNPS" selected="false">
            <parameter name="Meta data file" selected="false"/>
            <parameter name="Export ion identity networks">true</parameter>
            <parameter name="Presets">HIGHRES</parameter>
            <parameter name="Job title"/>
            <parameter name="Email"/>
            <parameter name="Username"/>
            <parameter name="Password"/>
            <parameter name="Open website">true</parameter>
        </parameter>
        <parameter name="Open folder">false</parameter>
    </batchstep>
</batch>
"""
