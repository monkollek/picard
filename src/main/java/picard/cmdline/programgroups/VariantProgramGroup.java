package picard.cmdline.programgroups;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;
import picard.util.help.HelpConstants;

/**
 * Program group for tools that manipulate variants and associated metadata
 */
public final class VariantProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return HelpConstants.DOC_CAT_VARIANT; }

    @Override
    public String getDescription() { return  HelpConstants.DOC_CAT_VARIANT_SUMMARY; }
}
