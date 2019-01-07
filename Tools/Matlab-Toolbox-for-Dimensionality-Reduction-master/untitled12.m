[mapped_data, mapping] = compute_mapping(NE_TSNE_Matrix(:,1:5), 'KernelPCA', 3);

scatter3(mapped_data(:,1),mapped_data(:,2),mapped_data(:,3),NE_Values_Top)