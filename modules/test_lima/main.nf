process TEST_LIMA {
    tag "Test Lima"
    publishDir "${params.outdir}/test_lima", mode: 'copy'
    
    output:
    path "lima_test.txt"

    script:
    """
    set -eu
    
    echo "=== LIMA TEST ===" > lima_test.txt
    echo "Apptainer version:" >> lima_test.txt
    apptainer --version >> lima_test.txt
    echo "" >> lima_test.txt
    echo "LIMA instances:" >> lima_test.txt
    limactl ls >> lima_test.txt
    echo "" >> lima_test.txt
    echo "=== TEST COMPLETE ===" >> lima_test.txt
    """

    stub:
    """
    echo "LIMA TEST - Apptainer version: 0.0.0" > lima_test.txt
    echo "LIMA instances:" >> lima_test.txt
    echo "NAME         STATUS     SSH                CPUS    MEMORY    DISK      DIR" >> lima_test.txt
    echo "apptainer    Running    123.456.789.0:49776    4       4GiB      100GiB    ~/.lima/apptainer" >> lima_test.txt
    """
}