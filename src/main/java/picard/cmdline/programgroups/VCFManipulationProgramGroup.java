package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Tools that manipulate variant call format (VCF) data
 */
public class VCFManipulationProgramGroup  implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_VCF_MANIPULATION; }

    @Override
    public String getDescription() { return HelpConstants.DOC_CAT_VCF_MANIPULATION_SUMMARY; }
}
