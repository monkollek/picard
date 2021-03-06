package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.*;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.EnumSet;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Simple little class that combines multiple VCFs that have exactly the same set of samples
 * and totally discrete sets of loci.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Gathers multiple VCF files from a scatter operation into a single VCF file. Input files " +
                "must be supplied in genomic order and must not have events at overlapping positions.",
        oneLineSummary = "Gathers multiple VCF files from a scatter operation into a single VCF file",
        programGroup = VcfOrBcf.class)
@DocumentedFeature
public class GatherVcfs extends CommandLineProgram {

    @Argument(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME,  doc="Input VCF file(s).")
	public List<File> INPUT;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF file.")
	public File OUTPUT;

    private static final Log log = Log.getInstance(GatherVcfs.class);

    public static void main(final String[] args) {
        new GatherVcfs().instanceMainWithExit(args);
    }

    public GatherVcfs() {
        CREATE_INDEX = true;
    }

    @Override
    protected int doWork() {
        log.info("Checking inputs.");
        INPUT = IOUtil.unrollFiles(INPUT, IOUtil.VCF_EXTENSIONS);
        for (final File f: INPUT) IOUtil.assertFileIsReadable(f);
        IOUtil.assertFileIsWritable(OUTPUT);

		final SAMSequenceDictionary sequenceDictionary = VCFFileReader.getSequenceDictionary(INPUT.get(0));

        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new PicardException("In order to index the resulting VCF input VCFs must contain ##contig lines.");
        }

        log.info("Checking file headers and first records to ensure compatibility.");
        assertSameSamplesAndValidOrdering(INPUT);

        if (areAllBlockCompressed(INPUT) && areAllBlockCompressed(CollectionUtil.makeList(OUTPUT))) {
            log.info("Gathering by copying gzip blocks. Will not be able to validate position non-overlap of files.");
            if (CREATE_INDEX) log.warn("Index creation not currently supported when gathering block compressed VCFs.");
            gatherWithBlockCopying(INPUT, OUTPUT);
        }
        else {
            log.info("Gathering by conventional means.");
            gatherConventionally(sequenceDictionary, CREATE_INDEX, INPUT, OUTPUT);
        }

        return 0;
    }

    /** Checks (via filename checking) that all files appear to be block compressed files. */
    private boolean areAllBlockCompressed(final List<File> input) {
        for (final File f : input) {
            if (VCFFileReader.isBCF(f) || !AbstractFeatureReader.hasBlockCompressedExtension(f)) return false;
        }

        return true;
    }

    /** Validates that all headers contain the same set of genotyped samples and that files are in order by position of first record. */
    private static void assertSameSamplesAndValidOrdering(final List<File> inputFiles) {
        final VCFHeader header = new VCFFileReader(inputFiles.get(0), false).getFileHeader();
        final SAMSequenceDictionary dict = header.getSequenceDictionary();
        final VariantContextComparator comparator = new VariantContextComparator(header.getSequenceDictionary());
        final List<String> samples = header.getGenotypeSamples();

        File lastFile = null;
        VariantContext lastContext = null;

        for (final File f : inputFiles) {
            final VCFFileReader in = new VCFFileReader(f, false);
            dict.assertSameDictionary(in.getFileHeader().getSequenceDictionary());
            final List<String> theseSamples = in.getFileHeader().getGenotypeSamples();

            if (!samples.equals(theseSamples)) {
                final SortedSet<String> s1 = new TreeSet<String>(samples);
                final SortedSet<String> s2 = new TreeSet<String>(theseSamples);
                s1.removeAll(theseSamples);
                s2.removeAll(samples);

                throw new IllegalArgumentException("VCFs do not have identical sample lists." +
                        " Samples unique to first file: " + s1 + ". Samples unique to " + f.getAbsolutePath() + ": " + s2 + ".");
            }

            final CloseableIterator<VariantContext> variantIterator = in.iterator();
            if (variantIterator.hasNext()) {
                final VariantContext currentContext = variantIterator.next();
                if (lastContext != null && comparator.compare(lastContext, currentContext) >= 0) {
                    throw new IllegalArgumentException("First record in file " + f.getAbsolutePath() + " is not after first record in " +
                            "previous file " + lastFile.getAbsolutePath());
                }

                lastContext = currentContext;
                lastFile    = f;
            }

            CloserUtil.close(in);
        }
    }

    /** Code for gathering multiple VCFs that works regardless of input format and output format, but can be slow. */
    private static void gatherConventionally(final SAMSequenceDictionary sequenceDictionary,
                                      final boolean createIndex,
                                      final List<File> inputFiles,
                                      final File outputFile) {
        final EnumSet<Options> options = EnumSet.copyOf(VariantContextWriterBuilder.DEFAULT_OPTIONS);
        if (createIndex){
            options.add(Options.INDEX_ON_THE_FLY);
        } else {
            options.remove(Options.INDEX_ON_THE_FLY);
        }
        final VariantContextWriter out = new VariantContextWriterBuilder()
                .setOptions(options)
                .setOutputFile(outputFile)
                .setReferenceDictionary(sequenceDictionary)
                .build();

        final ProgressLogger progress = new ProgressLogger(log, 10000);
        VariantContext lastContext = null;
        File lastFile = null;
        VCFHeader firstHeader = null;
        VariantContextComparator comparator = null;

        for (final File f : inputFiles) {
            log.debug("Gathering from file: ", f.getAbsolutePath());
            final VCFFileReader variantReader = new VCFFileReader(f, false);
            final PeekableIterator<VariantContext> variantIterator = new PeekableIterator<VariantContext>(variantReader.iterator());
            final VCFHeader header = variantReader.getFileHeader();

            if (firstHeader == null) {
                firstHeader = header;
                out.writeHeader(firstHeader);
                comparator = new VariantContextComparator(firstHeader.getContigLines());
            }

            if (lastContext != null && variantIterator.hasNext()) {
                final VariantContext vc = variantIterator.peek();
                if (comparator.compare(vc, lastContext) <= 0) {
                    throw new IllegalStateException("First variant in file " + f.getAbsolutePath() + " is at " + vc.getSource() +
                            " but last variant in earlier file " + lastFile.getAbsolutePath() + " is at " + lastContext.getSource());
                }
            }

            while (variantIterator.hasNext()) {
                lastContext = variantIterator.next();
                out.add(lastContext);
                progress.record(lastContext.getContig(), lastContext.getStart());
            }

            lastFile = f;

            CloserUtil.close(variantIterator);
            CloserUtil.close(variantReader);
        }

        out.close();
    }

    /**
     * Assumes that all inputs and outputs are block compressed VCF files and copies them without decompressing and parsing
     * most of the gzip blocks. Will decompress and parse blocks up to the one containing the end of the header in each file
     * (often the first block) and re-compress any data remaining in that block into a new block in the output file. Subsequent
     * blocks (excluding a terminator block if present) are copied directly from input to output.
     */
    private static void gatherWithBlockCopying(final List<File> vcfs, final File output) {
        try {
            final FileOutputStream out = new FileOutputStream(output);
            boolean isFirstFile = true;

            for (final File f : vcfs) {
                log.info("Gathering " + f.getAbsolutePath());
                final FileInputStream in = new FileInputStream(f);

                // a) It's good to check that the end of the file is valid and b) we need to know if there's a terminator block and not copy it
                final BlockCompressedInputStream.FileTermination term = BlockCompressedInputStream.checkTermination(f);
                if (term == BlockCompressedInputStream.FileTermination.DEFECTIVE) throw new PicardException(f.getAbsolutePath() + " does not have a valid GZIP block at the end of the file.");

                if (!isFirstFile) {
                    final BlockCompressedInputStream blockIn = new BlockCompressedInputStream(in, false);
                    boolean lastByteNewline = true;

                    while (blockIn.available() > 0) {
                        // Read a block - blockIn.available() is guaranteed to return the bytes remaining in the block that has been
                        // read, and since we haven't consumed any yet, that is the block size.
                        final int blockLength = blockIn.available();
                        final byte[] blockContents = new byte[blockLength];
                        final int read = blockIn.read(blockContents);
                        if (blockLength == 0 || read != blockLength) throw new IllegalStateException("Could not read available bytes from BlockCompressedInputStream.");

                        // Scan forward within the block to see if we can find the end of the header within this block
                        int firstNonHeaderByteIndex = -1;
                        for (int i=0; i<read; ++i) {
                            final byte b = blockContents[i];
                            final boolean thisByteNewline = (b == '\n' || b == '\r');

                            if (lastByteNewline && !thisByteNewline && b != '#') {
                                // Aha!  Found first byte of non-header data in file!
                                firstNonHeaderByteIndex = i;
                                break;
                            }

                            lastByteNewline = thisByteNewline;
                        }

                        // If we found the end of the header then write the remainder of this block out as a
                        // new gzip block and then break out of the while loop
                        if (firstNonHeaderByteIndex >= 0) {
                            final BlockCompressedOutputStream blockOut = new BlockCompressedOutputStream(out, null);
                            blockOut.write(blockContents, firstNonHeaderByteIndex, blockContents.length - firstNonHeaderByteIndex);
                            blockOut.flush();
                            // Don't close blockOut because closing underlying stream would break everything
                            break;
                        }
                    }
                }

                // Copy remainder of input stream into output stream
                final long currentPos = in.getChannel().position();
                final long length     = f.length();
                final long skipLast   = (term == BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK) ?
                        BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK.length : 0;
                final long bytesToWrite = length - skipLast - currentPos;

                IOUtil.transferByStream(in, out, bytesToWrite);
                in.close();
                isFirstFile = false;
            }

            // And lastly add the Terminator block and close up
            out.write(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
            out.close();
        }
        catch (final IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }
}
