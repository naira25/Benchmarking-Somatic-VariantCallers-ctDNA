BASE_DIR="/Users/nairaramosandres/Benchmarking-Somatic-VariantCallers-ctDNA"

subsamples=("1" "3" "5")
samples=("D" "E")

for sub in "${subsamples[@]}"; do
    VARIANTS_DIR="$BASE_DIR/Results/Results_Subsample_${sub}/variant_calling/benchmarking_vcf"
        
    cd "$VARIANTS_DIR" || continue

    for sam in "${samples[@]}"; do

        bcftools isec -n+2 -c all *_${sam}_*final.vcf.gz -p "ensemble_B_${sam}_2"
        for f in ensemble_B_${sam}_2/000*.vcf; do
            bgzip -f "$f"
            bcftools index -f "${f}.gz"
        done

        bcftools merge --force-samples -m none ensemble_B_${sam}_2/000*.vcf.gz | bcftools sort -Oz -o "S${sub}_B_${sam}_ensemble_2.vcf.gz"
        bcftools index -c "S${sub}_B_${sam}_ensemble_2.vcf.gz"

        bcftools isec -n=4 -c all *_${sam}_*final.vcf.gz -p "ensemble_B_${sam}_4"
        bgzip -c "ensemble_B_${sam}_4/0000.vcf" > "S${sub}_B_${sam}_ensemble_4.vcf.gz"
        bcftools index -c "S${sub}_B_${sam}_ensemble_4.vcf.gz"
        
    done

    cd "$BASE_DIR"
done