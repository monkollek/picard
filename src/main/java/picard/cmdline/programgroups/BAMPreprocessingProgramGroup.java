package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that align reads, flag duplicates and recalibrate base qualities
 */
public class BAMPreprocessingProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_BAM_PREPROCESSING; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_BAM_PREPROCESSING_SUMMARY; }
}